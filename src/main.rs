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

fn reduce_single_qual(q: u8, qual_reduction: &u8) -> u8 {
    let q_reduced: u8;
    if q >= *qual_reduction {
        q_reduced = q - qual_reduction;
    } else {
        q_reduced = 0;
    }
    q_reduced
}

fn reduce_qualities_in_edges(mut qual: Vec<u8>, num_bases: &usize, qual_reduction: &u8) -> Vec<u8> {
    let seq_len = qual.len() as usize;
    for pos in 0..*num_bases {
        if pos >= seq_len {
            continue;
        };
        qual[pos] = reduce_single_qual(qual[pos], qual_reduction);

        let pos_from_end = seq_len - pos - 1;
        if pos_from_end < *num_bases {
            continue;
        }
        qual[pos_from_end] = reduce_single_qual(qual[pos_from_end], qual_reduction);
    }
    qual
}

fn reduce_qualities_in_read(record: &mut Record, num_bases: &usize, qual_reduction: &u8) {
    let mut qual = record.qual().to_vec();

    qual = reduce_qualities_in_edges(qual, num_bases, qual_reduction);

    let cigar = record.cigar().take();
    let mut seq = record.seq().as_bytes();
    let qname = record.qname().to_owned();
    record.set(&qname, Some(&cigar), &mut seq, &qual);
}

fn trim_qualities_from_edges_in_bam(
    input_bam: &str,
    output_bam: &str,
    num_bases: &usize,
    qual_reduction: &u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = match input_bam {
        "-" => bam::Reader::from_stdin()?,
        _ => bam::Reader::from_path(input_bam)?,
    };

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

        reduce_qualities_in_read(&mut record, num_bases, qual_reduction);

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
    let num_bases: usize = 3;
    let qual_reduction: u8 = 20;

    // Call the quality reduction function
    trim_qualities_from_edges_in_bam(input_bam, output_bam, &num_bases, &qual_reduction)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_qualities_in_edges() {
        let qual = vec![30, 25, 20, 15, 10, 5];
        let num_bases = 2;
        let qual_reduction = 5;
        let expected = vec![25, 20, 20, 15, 5, 0];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reduce_qualities_in_edges_no_reduction() {
        let qual = vec![30, 25, 20, 15, 10, 5];
        let num_bases = 2;
        let qual_reduction = 0;
        let expected = vec![30, 25, 20, 15, 10, 5];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reduce_qualities_in_edges_full_reduction() {
        let qual = vec![30, 25, 20, 15, 10, 5];
        let num_bases = 3;
        let qual_reduction = 40;
        let expected = vec![0, 0, 0, 0, 0, 0];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reduce_qualities_in_overlaping_edges() {
        let qual = vec![30, 25, 20, 15, 10, 5];
        let num_bases = 4;
        let qual_reduction = 10;
        let expected = vec![20, 15, 10, 5, 0, 0];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reduce_qualities_in_edges_empty_qual() {
        let qual = vec![];
        let num_bases = 2;
        let qual_reduction = 5;
        let expected = vec![];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_reduce_qualities_in_edges_num_bases_exceeds_length() {
        let qual = vec![30, 25, 20];
        let num_bases = 4;
        let qual_reduction = 10;
        let expected = vec![20, 15, 10];
        let result = reduce_qualities_in_edges(qual, &num_bases, &qual_reduction);
        assert_eq!(result, expected);
    }
}
