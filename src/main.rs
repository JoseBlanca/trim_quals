use rust_htslib::{bam, bam::Format, bam::Header, bam::Read, bam::Record, bam::Writer};

fn reduce_quality(input_bam: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = bam::Reader::from_path(input_bam)?;
    //let mut reader = bam::Reader::from_stdin()?;

    let header_view = reader.header();
    let header = Header::from_template(header_view);

    // Open the output BAM file
    let format = Format::Cram;
    let format = Format::Sam;
    //let mut writer = bam::Writer::from_path(output_bam, &header)?;
    let mut writer = bam::Writer::from_stdout(&header, format)?;

    // Iterate through all records
    let mut record = Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.expect("Failed to parse record");
        // Write the modified record to the output BAM file
        writer.write(&record)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <input.bam>", args[0]);
        std::process::exit(1);
    }

    let input_bam = &args[1];

    // Call the quality reduction function
    reduce_quality(input_bam)?;

    Ok(())
}
