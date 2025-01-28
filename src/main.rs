use rust_htslib::{
    bam::{self, Format, Header, Read, Record},
    htslib::htsFormat,
};

fn get_file_format(hst_format: u32) -> Result<Format, String> {
    // formats (taken from htslib/hts.h enum htsExactFormat)
    match hst_format {
        3 => Ok(Format::Sam),
        4 => Ok(Format::Bam),
        6 => Ok(Format::Cram),
        _ => Err("Input file is not recognized as Sam, Bam or Cram".to_string()),
    }
}

fn reduce_qualities(record: &mut Record, num_bases: usize, qual_reduction: u8) {
    let mut qual = record.qual().to_vec();

    let seq_len = qual.len() as usize;
    for i in 0..num_bases {
        if i >= seq_len {
            continue;
        };
        let q = qual[i];
        let q_reduced: u8;
        if q >= qual_reduction {
            q_reduced = q - qual_reduction;
        } else {
            q_reduced = 0;
        }
        qual[i] = q_reduced;

        let i = i + 1;
        if i > seq_len {
            continue;
        }
        if seq_len - i <= i {
            continue;
        }
        let pos_from_end = seq_len - i;
        let q = qual[pos_from_end];
        let q_reduced: u8;
        if q >= qual_reduction {
            q_reduced = q - qual_reduction;
        } else {
            q_reduced = 0;
        }
        qual[pos_from_end] = q_reduced;
    }

    let cigar = record.cigar().take();
    let mut seq = record.seq().as_bytes();
    let qname = record.qname().to_owned();
    record.set(&qname, Some(&cigar), &mut seq, &qual);
}

fn reduce_quality(input_bam: &str, output_bam: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = match input_bam {
        "-" => bam::Reader::from_stdin()?,
        _ => bam::Reader::from_path(input_bam)?,
    };

    let num_bases: usize = 3;
    let qual_reduction: u8 = 20;

    let hst_format: htsFormat;
    unsafe {
        hst_format = (*reader.htsfile()).format;
    }
    // category 1 means sequence data
    if hst_format.category != 1 {
        return Err("The file is not recognized as sequence data".into());
    }
    let format = get_file_format(hst_format.format)?;

    let header_view = reader.header();
    let header = Header::from_template(header_view);

    let mut writer = match output_bam {
        "-" => bam::Writer::from_stdout(&header, format)?,
        _ => bam::Writer::from_path(output_bam, &header, format)?,
    };

    // Iterate through all records
    let mut record = Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.expect("Failed to parse record");

        reduce_qualities(&mut record, num_bases, qual_reduction);

        // Write the modified record to the output BAM file
        writer.write(&record)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.bam> <output.bam>", args[0]);
        std::process::exit(1);
    }

    let input_bam = &args[1];
    let output_bam = &args[2];

    // Call the quality reduction function
    reduce_quality(input_bam, output_bam)?;

    Ok(())
}
