extern crate clap;
extern crate env_logger;
extern crate exitcode;
extern crate log;

use clap::{Arg, App};
use log::{info, error};
use rust_htslib::{bam, bam::Read};

use std::collections::HashMap;

//all the relevant flags according to SAM specs and flagstat docs
//flagstat - http://www.htslib.org/doc/samtools-flagstat.html
const PAIRED_FLAG: u16 = 0x1;
const PROPERLY_ALIGNED_FLAG: u16 = 0x2;
const UNMAPPED_FLAG: u16 = 0x4;
const MATE_UNMAPPED_FLAG: u16 =  0x8;
const READ1_FLAG: u16 = 0x40;
const READ2_FLAG: u16 = 0x80;
const SECONDARY_FLAG: u16 = 0x100;
const FAIL_FLAG: u16 = 0x200;
const DUPLICATE_FLAG: u16 = 0x400;
const SUPPLEMENTARY_FLAG: u16 = 0x800;

//just rope this in for ease
const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");

fn main() {
    //initialize logging for our benefit later
    env_logger::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //open up our test file
    let matches = App::new("rust-flagstat")
        .version(VERSION.unwrap_or("?"))
        .author("J. Matthew Holt <jholt@hudsonalpha.org>")
        .about("flagstat - Rust implementation")
        .arg(Arg::with_name("in_fn")
            .short("i")
            .long("--input")
            .takes_value(true)
            .help("the BAM file to gather stats on"))
        .get_matches();
    
    let in_fn = match matches.value_of("in_fn") {
        Some(s) => s.to_string(),
        None => {
            error!("an input BAM file must be specified with --input");
            std::process::exit(exitcode::NOINPUT);
        }
    };
    let mut bam = match bam::Reader::from_path(&in_fn) {
        Ok(b) => b,
        Err(e) => {
            error!("Failed to open input file {:?}", in_fn);
            std::process::exit(exitcode::IOERR);
        }
    };

    let mut flag_counts = HashMap::new();
    
    let mut total_count: u64 = 0;
    let mut qc_failed: u64 = 0;
    let mut secondary_count: u64 = 0;
    let mut supplementary_count: u64 = 0;
    let mut duplicate_count: u64 = 0;
    let mut mapped_count: u64 = 0;
    let mut paired_count: u64 = 0;
    let mut read1_count: u64 = 0;
    let mut read2_count: u64 = 0;
    let mut properpair_count: u64 = 0;
    let mut bothmapped_count: u64 = 0;
    let mut singleton_count: u64 = 0;
    let mut different_chrom_count: u64 = 0;
    let mut hq_different_chrom_count: u64 = 0;
    
    //goes through and counts up each occurrence of tuple (flag-value, pair on same chrom, high quality mapping)
    for r in bam.records() {
        let record = r.unwrap();
        let flag = record.flags();
        let same_chrom = record.tid() == record.mtid();
        let hq_mapping = record.mapq() >= 5;

        let flag_entry = flag_counts.entry((flag, same_chrom, hq_mapping)).or_insert(0);
        *flag_entry += 1;
    }

    //now operate on the entries
    for ((flag, same_chrom, hq_mapping), fl_count) in &flag_counts {
        //total always goes up
        total_count += fl_count;

        //check if the qc is failed
        //TODO: if we care about breaking it into passed/failed like flagstat, then more needs to happen here
        if (flag & FAIL_FLAG) != 0 {
            qc_failed += fl_count;
        }

        //check if secondary alignment
        let is_secondary: bool = (flag & SECONDARY_FLAG) != 0;
        if is_secondary {
            secondary_count += fl_count;
        }

        //check if supplementary
        //NOTE: BUT it is apparently only added to counts if NOT secondary also
        let is_supplementary: bool = (flag & SUPPLEMENTARY_FLAG) != 0;
        if !is_secondary & is_supplementary {
            supplementary_count += fl_count;
        }

        //duplicate count
        //TODO: we didn't have any in our pipeline, may be more logic needed here
        if (flag & DUPLICATE_FLAG) != 0 {
            duplicate_count += fl_count;
        }

        //check for mapped
        let is_mapped: bool = (flag & UNMAPPED_FLAG) == 0;
        if is_mapped {
            mapped_count += fl_count;
        }

        //now a logic block for paired reads that are primary
        let is_paired: bool = (flag & PAIRED_FLAG) != 0;
        //TODO: i don't have any dups in here to test with, but we might need a is_duplicate check
        if is_paired & !is_secondary & !is_supplementary {
            //increment counts
            paired_count += fl_count;
            
            //count read1 and read2
            if (flag & READ1_FLAG) != 0 {
                read1_count += fl_count;
            }
            if (flag & READ2_FLAG) != 0 {
                read2_count += fl_count;
            }

            //logic block for primary, paired, mapped reads
            if is_mapped {
                //properly paired counts
                if (flag & PROPERLY_ALIGNED_FLAG) != 0 {
                    properpair_count += fl_count;
                }

                //checks when both this read and the pair are mapped
                let is_matemapped: bool = (flag & MATE_UNMAPPED_FLAG) == 0;
                if is_matemapped {
                    bothmapped_count += fl_count;

                    //different chromosome mappings and quality of mapping comes into play here only
                    if !same_chrom {
                        different_chrom_count += fl_count;
                        if *hq_mapping {
                            hq_different_chrom_count += fl_count;
                        }
                    }
                } else {
                    singleton_count += fl_count;
                }
            }
        }
    }
    
    //print out something similar to flagstat, but really just want to check numbers are same
    println!("total_count: {:?}", total_count);
    println!("qc_failed: {:?}", qc_failed);
    println!("secondary_count: {:?}", secondary_count);
    println!("supplementary_count: {:?}", supplementary_count);
    println!("duplicate_count: {:?}", duplicate_count);
    println!("mapped_count: {:?}", mapped_count);
    println!("paired_count: {:?}", paired_count);
    println!("read1: {:?}", read1_count);
    println!("read2: {:?}", read2_count);
    println!("properpair_count: {:?}", properpair_count);
    println!("bothmapped_count: {:?}", bothmapped_count);
    println!("singleton_count: {:?}", singleton_count);
    println!("different_chrom_count: {:?}", different_chrom_count);
    println!("hq_different_chrom_count: {:?}", hq_different_chrom_count);
}
