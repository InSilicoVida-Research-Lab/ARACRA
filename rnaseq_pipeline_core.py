#!/usr/bin/env python3
"""
RNA-seq Processing Pipeline - Core Module
Core functionality for processing RNA sequencing data from SRA
"""

import os
import subprocess
import pandas as pd
import threading
import time
from pathlib import Path
import shutil
import glob
import traceback
import json
import logging
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import psutil
import tempfile
import re
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import queue

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class PipelineConfig:
    """Configuration for pipeline parameters"""
    output_dir: str
    trimmed_dir: str
    hisat2_dir: str
    trimmomatic_jar: str
    adapters: str
    index_path: str
    reference_gtf: str
    num_processes: int = 4
    featurecounts_threads: int = 8
    max_retries: int = 3
    timeout_prefetch: int = 1800
    timeout_fasterq: int = 3600
    timeout_trim: int = 1800
    timeout_align: int = 3600
    timeout_count: int = 1800
    cleanup_intermediates: bool = True
    
    def validate(self) -> Tuple[bool, List[str]]:
        """Validate configuration parameters"""
        errors = []
        
        # Check required files
        files_to_check = {
            'Trimmomatic JAR': self.trimmomatic_jar,
            'Adapter sequences': self.adapters,
            'Reference GTF': self.reference_gtf
        }
        
        for name, path in files_to_check.items():
            if not path or not os.path.exists(path):
                errors.append(f"{name} not found: {path}")
        
        # Validate HISAT2 index
        if not self.index_path:
            errors.append("HISAT2 index path not specified")
        else:
            is_valid, message = validate_hisat2_index(self.index_path)
            if not is_valid:
                errors.append(f"HISAT2 index error: {message}")
        
        # Check numeric parameters
        cpu_count = os.cpu_count() or 1
        if self.num_processes < 1 or self.num_processes > cpu_count:
            errors.append(f"Invalid number of processes: {self.num_processes}")
        
        return len(errors) == 0, errors

@dataclass
class ProcessingResult:
    """Result for individual accession processing"""
    accession: str
    status: str  # 'completed', 'failed', 'timeout', 'cancelled'
    start_time: float
    end_time: Optional[float] = None
    error: Optional[str] = None
    steps: Dict[str, str] = None
    processing_time: Optional[float] = None
    retry_count: int = 0
    
    def __post_init__(self):
        if self.steps is None:
            self.steps = {}
        if self.end_time and self.start_time:
            self.processing_time = self.end_time - self.start_time

class ThreadSafeLogger:
    """Thread-safe logger for pipeline messages"""
    
    def __init__(self, max_messages: int = 1000):
        self._lock = threading.RLock()
        self.max_messages = max_messages
        self._messages = []
        
    def log_message(self, message: str, level: str = "INFO"):
        """Thread-safe logging"""
        try:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            formatted_message = f"[{timestamp}] [{level}] {message}"
            
            with self._lock:
                self._messages.append(formatted_message)
                if len(self._messages) > self.max_messages:
                    self._messages = self._messages[-self.max_messages:]
                
                # Log to standard logging
                getattr(logger, level.lower(), logger.info)(message)
                
        except Exception as e:
            print(f"Logging error: {e}")
    
    def get_messages(self) -> List[str]:
        """Get all log messages"""
        with self._lock:
            return self._messages.copy()

class ResourceMonitor:
    """Monitor system resources during processing"""
    
    def __init__(self):
        self.start_time = time.time()
        self.peak_memory = 0
        self.peak_cpu = 0
        self._monitoring = False
        self._monitor_thread = None
        self._stats = {
            'current_memory': 0,
            'current_cpu': 0,
            'peak_memory': 0,
            'peak_cpu': 0
        }
    
    def start_monitoring(self):
        """Start resource monitoring"""
        self._monitoring = True
        self._monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._monitor_thread.start()
    
    def stop_monitoring(self):
        """Stop resource monitoring"""
        self._monitoring = False
        if self._monitor_thread and self._monitor_thread.is_alive():
            self._monitor_thread.join(timeout=1)
    
    def _monitor_loop(self):
        """Monitor system resources"""
        while self._monitoring:
            try:
                process = psutil.Process()
                memory_mb = process.memory_info().rss / 1024 / 1024
                cpu_percent = process.cpu_percent()
                
                self._stats['current_memory'] = memory_mb
                self._stats['current_cpu'] = cpu_percent
                self._stats['peak_memory'] = max(self._stats['peak_memory'], memory_mb)
                self._stats['peak_cpu'] = max(self._stats['peak_cpu'], cpu_percent)
                
                time.sleep(5)
            except Exception:
                time.sleep(5)
    
    def get_stats(self) -> Dict[str, float]:
        """Get current resource statistics"""
        return self._stats.copy()

def validate_hisat2_index(index_path: str) -> Tuple[bool, str]:
    """Validate HISAT2 index files"""
    try:
        if not index_path or index_path.strip() == "":
            return False, "Index path not specified"
        
        # Handle different path formats
        if os.path.isfile(index_path):
            index_dir = os.path.dirname(index_path)
            index_base = os.path.splitext(os.path.basename(index_path))[0]
            full_base = index_path
        elif os.path.isdir(index_path):
            ht2_files = glob.glob(os.path.join(index_path, "*.ht2"))
            if ht2_files:
                first_file = os.path.basename(ht2_files[0])
                index_base = first_file.split('.')[0]
                full_base = os.path.join(index_path, index_base)
            else:
                return False, f"No .ht2 files found in directory: {index_path}"
        else:
            index_dir = os.path.dirname(index_path)
            index_base = os.path.basename(index_path)
            full_base = index_path
            
            if not index_dir:
                index_dir = os.getcwd()
                full_base = os.path.join(index_dir, index_base)
        
        # Check for required index files
        required_extensions = [f"{i}.ht2" for i in range(1, 9)]
        found_files = []
        
        for ext in required_extensions:
            index_file = f"{full_base}.{ext}"
            if os.path.exists(index_file):
                found_files.append(ext)
        
        if len(found_files) >= 8:
            return True, f"Found complete index with {len(found_files)} files"
        else:
            return False, f"Incomplete index - found {len(found_files)} files, expected 8"
            
    except Exception as e:
        return False, f"Validation error: {str(e)}"

def validate_tools() -> List[str]:
    """Check for required tools"""
    tools_config = {
        'prefetch': {'cmd': ['--version'], 'timeout': 10},
        'fasterq-dump': {'cmd': ['--version'], 'timeout': 10},
        'hisat2': {'cmd': ['--version'], 'timeout': 10},
        'samtools': {'cmd': ['--version'], 'timeout': 10},
        'featureCounts': {'cmd': ['-v'], 'timeout': 10},
        'java': {'cmd': ['-version'], 'timeout': 5}
    }
    
    missing_tools = []
    
    for tool, config in tools_config.items():
        try:
            result = subprocess.run(
                [tool] + config['cmd'],
                capture_output=True,
                text=True,
                timeout=config['timeout'],
                errors='replace'
            )
            
            if result.returncode != 0 and result.returncode != 1:
                if result.returncode == 127:  # Command not found
                    missing_tools.append(tool)
                    
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            missing_tools.append(tool)
        except Exception as e:
            logger.warning(f"Error checking tool {tool}: {e}")
            missing_tools.append(tool)
    
    return missing_tools

@contextmanager
def temp_directory():
    """Context manager for temporary directory"""
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def create_directories(directories: List[str]) -> Tuple[List[str], List[str]]:
    """Create required directories"""
    created_dirs = []
    failed_dirs = []
    
    for directory in directories:
        if not directory or directory.strip() == "":
            continue
        
        try:
            dir_path = Path(directory)
            dir_path.mkdir(parents=True, exist_ok=True)
            
            # Verify directory is writable
            test_file = dir_path / ".write_test"
            try:
                test_file.touch()
                test_file.unlink()
                created_dirs.append(str(dir_path))
            except PermissionError:
                failed_dirs.append(f"{directory} (not writable)")
                
        except Exception as e:
            failed_dirs.append(f"{directory} ({str(e)})")
    
    return created_dirs, failed_dirs

def run_command_with_retry(cmd: List[str], timeout: int, max_retries: int = 3, 
                          cwd: Optional[str] = None) -> subprocess.CompletedProcess:
    """Run command with retry logic"""
    for attempt in range(max_retries):
        try:
            logger.info(f"Running command (attempt {attempt + 1}/{max_retries}): {' '.join(cmd[:3])}...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=cwd,
                errors='replace'
            )
            
            if result.returncode == 0:
                return result
            else:
                logger.warning(f"Command failed with code {result.returncode}")
                if attempt == max_retries - 1:
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
                    
        except subprocess.TimeoutExpired:
            logger.warning(f"Command timeout (attempt {attempt + 1}/{max_retries})")
            if attempt == max_retries - 1:
                raise
        except Exception as e:
            logger.warning(f"Command error (attempt {attempt + 1}/{max_retries}): {e}")
            if attempt == max_retries - 1:
                raise
        
        time.sleep(2 ** attempt)  # Exponential backoff
    
    raise Exception("Max retries exceeded")

def run_combined_featurecounts(config: PipelineConfig, bam_files: List[str], 
                              thread_logger: ThreadSafeLogger) -> bool:
    """Run featureCounts on all BAM files together to create a combined count matrix"""
    try:
        if not bam_files:
            thread_logger.log_message("No BAM files to process for feature counting", "WARNING")
            return False
        
        thread_logger.log_message(f"Running combined feature counting for {len(bam_files)} samples")
        
        # Output file for combined counts - using .out format as requested
        combined_counts_file = os.path.join(config.hisat2_dir, "combined_counts.out")
        
        # Build featureCounts command for single-end data
        featurecounts_cmd = [
            "featureCounts",
            "-a", config.reference_gtf,
            "-o", combined_counts_file,
            "-T", str(config.featurecounts_threads),
            "-g", "gene_id",
            "-t", "exon",
            # REMOVED -p, -B, -C flags (these are for paired-end only)
            # Added essential flags for single-end RNA-seq
            "-s", "0",  # Unstranded library
            "-Q", "10",  # Minimum mapping quality
            "-M",  # Count multi-mapping reads
            "--primary",  # Count primary alignments only
        ] + bam_files
        
        thread_logger.log_message(f"Running featureCounts with {len(bam_files)} BAM files (single-end mode)")
        
        # Run featureCounts
        result = subprocess.run(
            featurecounts_cmd,
            capture_output=True,
            text=True,
            timeout=config.timeout_count * len(bam_files),
            errors='replace'
        )
        
        if result.returncode == 0:
            thread_logger.log_message(f"Combined feature counting completed successfully")
            thread_logger.log_message(f"Count matrix saved to: {combined_counts_file}")
            
            # Log the summary statistics
            summary_file = combined_counts_file + ".summary"
            if os.path.exists(summary_file):
                try:
                    with open(summary_file, 'r') as f:
                        summary_lines = f.readlines()[:20]  # First 20 lines should show assignment stats
                        thread_logger.log_message("FeatureCounts assignment summary:")
                        for line in summary_lines:
                            thread_logger.log_message(line.strip())
                except Exception as e:
                    thread_logger.log_message(f"Could not read summary file: {e}", "WARNING")
            
            return True
        else:
            thread_logger.log_message(f"FeatureCounts failed with return code {result.returncode}", "ERROR")
            thread_logger.log_message(f"Error output: {result.stderr}", "ERROR")
            return False
            
    except subprocess.TimeoutExpired:
        thread_logger.log_message("FeatureCounts timeout exceeded", "ERROR")
        return False
    except Exception as e:
        thread_logger.log_message(f"Error during combined feature counting: {str(e)}", "ERROR")
        return False

def process_single_accession(accession: str, config: PipelineConfig, 
                           thread_logger: ThreadSafeLogger) -> ProcessingResult:
    """Process a single SRA accession"""
    result = ProcessingResult(
        accession=accession,
        status="started",
        start_time=time.time()
    )
    
    try:
        thread_logger.log_message(f"Starting processing for {accession}")
        
        # Check if processing should be cancelled
        if os.path.exists(os.path.join(config.output_dir, ".stop_pipeline")):
            result.status = "cancelled"
            result.end_time = time.time()
            return result
        
        # Step 1: Download with prefetch
        thread_logger.log_message(f"Step 1/4: Downloading {accession}")
        
        prefetch_cmd = [
            "prefetch", accession, 
            "--output-directory", config.output_dir,
            "--max-size", "50G",
            "--progress"
        ]
        
        try:
            run_command_with_retry(prefetch_cmd, config.timeout_prefetch, config.max_retries, config.output_dir)
            result.steps["download"] = "success"
            thread_logger.log_message(f"Download completed for {accession}")
        except Exception as e:
            # Check if file already exists
            possible_locations = [
                os.path.join(config.output_dir, accession, f"{accession}.sra"),
                os.path.join(config.output_dir, f"{accession}.sra")
            ]
            
            if any(os.path.exists(loc) for loc in possible_locations):
                thread_logger.log_message(f"SRA file already exists for {accession}, continuing...")
                result.steps["download"] = "success"
            else:
                raise Exception(f"Prefetch failed: {e}")
        
        # Find SRA file
        sra_file = None
        for loc in [
            os.path.join(config.output_dir, accession, f"{accession}.sra"),
            os.path.join(config.output_dir, f"{accession}.sra")
        ]:
            if os.path.exists(loc):
                sra_file = loc
                break
        
        # Step 2: Convert to FASTQ
        thread_logger.log_message(f"Step 2/4: Converting {accession} to FASTQ")
        
        # Create a unique subdirectory for this accession's FASTQ output
        accession_fastq_dir = os.path.join(config.output_dir, f"{accession}_fastq_temp")
        os.makedirs(accession_fastq_dir, exist_ok=True)
        
        fasterq_cmd = [
            "fasterq-dump",
            sra_file if sra_file else accession,
            "--outdir", accession_fastq_dir,  # Use unique directory
            "--threads", "4",
            "--split-files",
            "--progress"
        ]
        
        with temp_directory() as temp_dir:
            fasterq_cmd.extend(["--temp", temp_dir])
            run_command_with_retry(fasterq_cmd, config.timeout_fasterq, config.max_retries)
        
        # Move FASTQ files to main output directory
        for filename in os.listdir(accession_fastq_dir):
            if filename.endswith(".fastq"):
                src = os.path.join(accession_fastq_dir, filename)
                dst = os.path.join(config.output_dir, filename)
                shutil.move(src, dst)
        
        # Remove temporary directory
        shutil.rmtree(accession_fastq_dir, ignore_errors=True)
        
        result.steps["convert"] = "success"
        thread_logger.log_message(f"FASTQ conversion completed for {accession}")
        
        # Find FASTQ files
        fastq_patterns = [
            os.path.join(config.output_dir, f"{accession}.fastq"),
            os.path.join(config.output_dir, f"{accession}_1.fastq"),
            os.path.join(config.output_dir, f"{accession}_2.fastq")
        ]
        
        input_fastq_files = [p for p in fastq_patterns if os.path.exists(p)]
        
        if not input_fastq_files:
            raise Exception(f"No FASTQ files found for {accession}")
        
        # Determine if paired-end
        is_paired = len(input_fastq_files) >= 2
        
        # Step 3: Trimming
        thread_logger.log_message(f"Step 3/4: Trimming {accession} ({'paired' if is_paired else 'single'}-end)")
        
        if is_paired:
            input_r1 = os.path.join(config.output_dir, f"{accession}_1.fastq")
            input_r2 = os.path.join(config.output_dir, f"{accession}_2.fastq")
            output_r1_paired = os.path.join(config.trimmed_dir, f"{accession}_1_paired.fastq")
            output_r1_unpaired = os.path.join(config.trimmed_dir, f"{accession}_1_unpaired.fastq")
            output_r2_paired = os.path.join(config.trimmed_dir, f"{accession}_2_paired.fastq")
            output_r2_unpaired = os.path.join(config.trimmed_dir, f"{accession}_2_unpaired.fastq")
            
            trim_cmd = [
                "java", "-Xmx4G", "-jar", config.trimmomatic_jar, "PE", 
                "-threads", "4", "-phred33",
                input_r1, input_r2,
                output_r1_paired, output_r1_unpaired,
                output_r2_paired, output_r2_unpaired,
                f"ILLUMINACLIP:{config.adapters}:2:30:10",
                "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
        else:
            input_fastq = input_fastq_files[0]
            output_fastq = os.path.join(config.trimmed_dir, f"{accession}_trimmed.fastq")
            
            trim_cmd = [
                "java", "-Xmx4G", "-jar", config.trimmomatic_jar, "SE",
                "-threads", "4", "-phred33",
                input_fastq, output_fastq,
                f"ILLUMINACLIP:{config.adapters}:2:30:10",
                "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
        
        run_command_with_retry(trim_cmd, config.timeout_trim, config.max_retries)
        result.steps["trim"] = "success"
        thread_logger.log_message(f"Trimming completed for {accession}")
        
        # Step 4: Alignment
        thread_logger.log_message(f"Step 4/4: Aligning {accession}")
        
        sam_file = os.path.join(config.hisat2_dir, f"{accession}.sam")
        
        if is_paired:
            hisat2_cmd = [
                "hisat2", "-x", config.index_path,
                "-1", os.path.join(config.trimmed_dir, f"{accession}_1_paired.fastq"),
                "-2", os.path.join(config.trimmed_dir, f"{accession}_2_paired.fastq"),
                "-S", sam_file, "--threads", "4", "--summary-file", 
                os.path.join(config.hisat2_dir, f"{accession}_summary.txt")
            ]
        else:
            hisat2_cmd = [
                "hisat2", "-x", config.index_path,
                "-U", os.path.join(config.trimmed_dir, f"{accession}_trimmed.fastq"),
                "-S", sam_file, 
                "--threads", "4",
                "--dta",  # Important for RNA-seq
                "--summary-file", os.path.join(config.hisat2_dir, f"{accession}_summary.txt")
            ]
        
        # Run HISAT2 alignment
        hisat2_result = run_command_with_retry(hisat2_cmd, config.timeout_align, config.max_retries)
        
        # Check if SAM file was created and has content
        if not os.path.exists(sam_file) or os.path.getsize(sam_file) == 0:
            raise Exception(f"HISAT2 alignment failed - SAM file not created or empty for {accession}")
        
        # Convert to sorted BAM
        sorted_bam = os.path.join(config.hisat2_dir, f"{accession}_sorted.bam")
        
        thread_logger.log_message(f"Converting SAM to sorted BAM for {accession}")
        
        sam_to_bam_cmd = [
            "samtools", "sort", 
            "-@", "4",
            "-o", sorted_bam,
            sam_file
        ]
        
        try:
            subprocess.run(sam_to_bam_cmd, check=True, timeout=config.timeout_align, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            thread_logger.log_message(f"samtools sort error: {e.stderr}", "ERROR")
            raise
        
        # Verify BAM file was created
        if not os.path.exists(sorted_bam) or os.path.getsize(sorted_bam) == 0:
            raise Exception(f"BAM file creation failed for {accession}")
        
        # Index BAM
        thread_logger.log_message(f"Indexing BAM file for {accession}")
        index_cmd = ["samtools", "index", sorted_bam]
        subprocess.run(index_cmd, check=True, timeout=600)
        
        result.steps["align"] = "success"
        thread_logger.log_message(f"Alignment completed for {accession}")
        
        # Cleanup intermediate files if requested
        if config.cleanup_intermediates:
            # Only cleanup files, not directories
            cleanup_patterns = [
                sam_file,  # SAM file
                os.path.join(config.output_dir, f"{accession}.fastq"),
                os.path.join(config.output_dir, f"{accession}_1.fastq"),
                os.path.join(config.output_dir, f"{accession}_2.fastq"),
            ]
            
            for file_path in cleanup_patterns:
                try:
                    if os.path.exists(file_path) and os.path.isfile(file_path):
                        os.remove(file_path)
                        thread_logger.log_message(f"Removed intermediate file: {file_path}")
                except Exception as e:
                    thread_logger.log_message(f"Warning: Could not remove {file_path}: {e}", "WARNING")
            
            # Optionally remove SRA directory
            sra_dir = os.path.join(config.output_dir, accession)
            try:
                if os.path.exists(sra_dir) and os.path.isdir(sra_dir):
                    shutil.rmtree(sra_dir)
                    thread_logger.log_message(f"Removed SRA directory: {sra_dir}")
            except Exception as e:
                thread_logger.log_message(f"Warning: Could not remove directory {sra_dir}: {e}", "WARNING")
        
        result.status = "completed"
        result.end_time = time.time()
        
        if result.processing_time:
            thread_logger.log_message(f"Successfully completed {accession} in {result.processing_time:.1f}s")
        else:
            thread_logger.log_message(f"Successfully completed {accession}")
        
        return result
        
    except Exception as e:
        result.status = "failed"
        result.error = str(e)
        result.end_time = time.time()
        thread_logger.log_message(f"Processing failed for {accession}: {e}", "ERROR")
        return result

def run_pipeline_thread(accession_numbers: List[str], config: PipelineConfig):
    """Main pipeline execution thread with parallel processing"""
    thread_logger = ThreadSafeLogger()
    
    try:
        thread_logger.log_message("Pipeline thread started")
        
        # Start resource monitoring
        resource_monitor = ResourceMonitor()
        resource_monitor.start_monitoring()
        
        successful = 0
        failed = 0
        detailed_results = []
        successful_bam_files = []  # Track successfully created BAM files
        total_accessions = len(accession_numbers)
        start_time = time.time()
        
        thread_logger.log_message(f"Processing {total_accessions} accessions with {config.num_processes} parallel processes")
        
        # Process accessions in parallel
        with ThreadPoolExecutor(max_workers=config.num_processes) as executor:
            # Submit all tasks
            future_to_accession = {
                executor.submit(process_single_accession, acc, config, thread_logger): acc 
                for acc in accession_numbers
            }
            
            # Process results as they complete
            for i, future in enumerate(as_completed(future_to_accession)):
                accession = future_to_accession[future]
                
                try:
                    # Check for cancellation
                    if os.path.exists(os.path.join(config.output_dir, ".stop_pipeline")):
                        thread_logger.log_message("Pipeline stopped by user")
                        executor.shutdown(wait=False)
                        break
                    
                    result = future.result()
                    detailed_results.append(result)
                    
                    if result.status == "completed":
                        successful += 1
                        # Track the BAM file for this successful sample
                        bam_file = os.path.join(config.hisat2_dir, f"{accession}_sorted.bam")
                        if os.path.exists(bam_file):
                            successful_bam_files.append(bam_file)
                            thread_logger.log_message(f"Completed {accession} successfully ({successful}/{total_accessions})")
                    else:
                        failed += 1
                        thread_logger.log_message(f"Failed processing {accession} ({failed} failed so far)", "WARNING")
                    
                    # Calculate progress
                    progress = (i + 1) / total_accessions
                    elapsed_time = time.time() - start_time
                    
                    if progress > 0:
                        estimated_total_time = elapsed_time / progress
                        estimated_remaining = estimated_total_time - elapsed_time
                    else:
                        estimated_remaining = 0
                    
                    # Store results
                    results_summary = {
                        'total_processed': i + 1,
                        'successful': successful,
                        'failed': failed,
                        'detailed_results': detailed_results,
                        'progress': progress,
                        'estimated_remaining': estimated_remaining,
                        'current_accession': f"Processing {config.num_processes} in parallel",
                        'processing_status': 'running'
                    }
                    
                    # Write results to file for main thread
                    try:
                        results_file = os.path.join(config.output_dir, ".pipeline_results.json")
                        with open(results_file, 'w') as f:
                            json.dump(results_summary, f, default=str)
                    except:
                        pass
                    
                except Exception as e:
                    error_msg = f"Critical error processing {accession}: {str(e)}"
                    thread_logger.log_message(error_msg, "ERROR")
                    failed += 1
                    
                    # Create a failed result
                    failed_result = ProcessingResult(
                        accession=accession,
                        status="failed",
                        start_time=time.time(),
                        end_time=time.time(),
                        error=str(e)
                    )
                    detailed_results.append(failed_result)
        
        thread_logger.log_message(f"All parallel processing completed. Successful: {successful}, Failed: {failed}")
        
        # Run combined feature counting on all successful BAM files
        if successful_bam_files:
            thread_logger.log_message(f"Starting combined feature counting for {len(successful_bam_files)} samples")
            featurecounts_success = run_combined_featurecounts(config, successful_bam_files, thread_logger)
            
            if featurecounts_success:
                thread_logger.log_message("Combined feature counting completed successfully")
            else:
                thread_logger.log_message("Combined feature counting failed", "ERROR")
        else:
            thread_logger.log_message("No successful alignments to process for feature counting", "WARNING")
        
        # Pipeline completion
        end_time = time.time()
        total_time = end_time - start_time
        
        final_results = {
            'total_processed': len(detailed_results),
            'successful': successful,
            'failed': failed,
            'detailed_results': detailed_results,
            'total_time': total_time,
            'completion_time': datetime.now().isoformat(),
            'processing_status': 'completed',
            'bam_files_processed': len(successful_bam_files),
            'count_matrix_generated': len(successful_bam_files) > 0
        }
        
        # Write final results
        try:
            results_file = os.path.join(config.output_dir, ".pipeline_results.json")
            with open(results_file, 'w') as f:
                json.dump(final_results, f, default=str)
        except:
            pass
        
        # Stop resource monitoring
        resource_monitor.stop_monitoring()
        
        thread_logger.log_message(f"Pipeline completed: {successful} successful, {failed} failed in {total_time:.1f}s")
        if len(successful_bam_files) > 0:
            thread_logger.log_message(f"Generated combined count matrix for {len(successful_bam_files)} samples")
        
    except Exception as e:
        error_msg = f"Critical pipeline error: {str(e)}\n{traceback.format_exc()}"
        thread_logger.log_message(error_msg, "ERROR")
        
        if 'resource_monitor' in locals():
            resource_monitor.stop_monitoring()

def export_results_to_csv(results: Dict) -> str:
    """Export detailed results to CSV format"""
    if not results.get('detailed_results'):
        return ""
    
    # Convert results to DataFrame
    rows = []
    for result in results['detailed_results']:
        if isinstance(result, ProcessingResult):
            row = asdict(result)
        elif isinstance(result, dict):
            row = result.copy()
        else:
            # Skip non-dict, non-ProcessingResult items
            continue
        
        # Flatten steps dictionary if it exists
        if 'steps' in row and isinstance(row['steps'], dict):
            steps = row.pop('steps', {})
            for step_name, step_status in steps.items():
                row[f'step_{step_name}'] = step_status
        
        rows.append(row)
    
    if not rows:
        return ""
    
    df = pd.DataFrame(rows)
    return df.to_csv(index=False)

def create_processing_summary(results: Dict) -> str:
    """Create a text summary of processing results"""
    if not results:
        return "No results available"
    
    summary = []
    summary.append("RNA-seq Processing Pipeline Summary")
    summary.append("=" * 40)
    summary.append(f"Total accessions processed: {results.get('total_processed', 0)}")
    summary.append(f"Successful: {results.get('successful', 0)}")
    summary.append(f"Failed: {results.get('failed', 0)}")
    
    if 'total_time' in results:
        summary.append(f"Total processing time: {results['total_time']:.1f} seconds")
    
    if 'completion_time' in results:
        summary.append(f"Completed at: {results['completion_time']}")
    
    # Add failed accessions details
    failed_accessions = []
    for result in results.get('detailed_results', []):
        if isinstance(result, dict) and result.get('status') == 'failed':
            failed_accessions.append(f"  • {result['accession']}: {result.get('error', 'Unknown error')}")
        elif hasattr(result, 'status') and result.status == 'failed':
            failed_accessions.append(f"  • {result.accession}: {result.error or 'Unknown error'}")
    
    if failed_accessions:
        summary.append("\nFailed accessions:")
        summary.extend(failed_accessions)
    
    return "\n".join(summary)

def validate_accession_format(accessions: List[str]) -> Tuple[List[str], List[str]]:
    """Validate SRA accession number format"""
    valid_accessions = []
    invalid_accessions = []
    
    sra_pattern = re.compile(r'^[SED]RR\d+')
    for acc in accessions:
        acc = acc.strip()
        if acc and sra_pattern.match(acc):
            valid_accessions.append(acc)
        elif acc:
            invalid_accessions.append(acc)
    
    return valid_accessions, invalid_accessions
