require 'fileutils'
require 'spec_helper'
require 'ribopip/array_writer'
require 'ribopip/counts'

counts = Ribopip::Counts::Bookkeeper.new
counts.parse_nreads('./spec/testdata/SRR315601_625.fastq', 'n_total')

describe Ribopip::ArrayWriter::TSVWriter do
  it '.write' do
    tsv = Ribopip::ArrayWriter::TSVWriter.new(counts.array, true)
    tsv.write('./spec/testdata/array_writer/SRR315601_625.counts.tsv')
    expect(
      FileUtils.compare_file(
        './spec/testdata/array_writer/SRR315601_625.counts.tsv',
        './spec/testdata/array_writer/expected_read_counts.tsv'
      )
    ).to be_truthy
  end
end

describe Ribopip::ArrayWriter::XMLWriter do
  it '.write' do
    xml = Ribopip::ArrayWriter::XMLWriter.new(counts.array, true)
    xml.write('./spec/testdata/array_writer/SRR315601_625.counts.xml')
    expect(
      FileUtils.compare_file(
        './spec/testdata/array_writer/SRR315601_625.counts.xml',
        './spec/testdata/array_writer/expected_read_counts.xml'
      )
    ).to be_truthy
  end
end

describe Ribopip::ArrayWriter::Writer do
  it '.write' do
    xml = Ribopip::ArrayWriter::Writer.new(counts.array, true)
    xml.write('./spec/testdata/array_writer/SRR315601_625.counts.txt')
    expect(
      FileUtils.compare_file(
        './spec/testdata/array_writer/SRR315601_625.counts.txt',
        './spec/testdata/array_writer/expected_read_counts.txt'
      )
    ).to be_truthy
  end
end
