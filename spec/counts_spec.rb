require 'spec_helper'
require 'ribopip/counts'

include Ribopip::Counts

describe Ribopip::Counts::Array do
  let(:arr) { Ribopip::Counts::Array.new }
  it '.new' do
    expect(arr.array).to be_empty
  end

  it '.push' do
    arr.push('n_total', 625)
    expect(arr.array).to eq [['n_total', 625, 100]]
    arr.push('n_foo', 125)
    expect(arr.array).to eq [['n_total', 625, 100], ['n_foo', 125, 20]]
  end

  it '.each' do
    yielded = []
    arr.push('n_total', 625)
    arr.push('n_foo', 125)
    arr.push('n_bar', 25)
    arr.each { |x| yielded << x }
    expect(yielded).to eq [['n_total', 625, 100],
                           ['n_foo', 125, 20],
                           ['n_bar', 25, 4]]
  end
end

describe Ribopip::Counts do
  it 'FLAGSTAT' do
    num = FLAGSTAT.call('./spec/testdata/counts/uni_mapped_0err.bam', nil)
    expect(num).to eq(166)
  end

  it 'LINECOUNT' do
    num = LINECOUNT.call('./spec/testdata/SRR315601_625.fastq', nil)
    expect(num).to eq(625)
  end

  it 'REGEX' do
    num = REGEX.call('./spec/testdata/counts/align_summary.txt', 'Mapped')
    expect(num).to eq(313)
  end
end

describe Ribopip::Counts::Bookkeeper do
  let(:counts) { Ribopip::Counts::Bookkeeper.new }

  it '.new' do
    yielded = []
    counts.array.each { |x| yielded << x }
    expect(yielded).to be_empty
  end

  it '.parse_nreads' do
    counts.parse_nreads('./spec/testdata/SRR315601_625.fastq',
                        'n_total')
    counts.parse_nreads('./spec/testdata/counts/align_summary.txt',
                        'n_mapped', 'Mapped')
    counts.parse_nreads('./spec/testdata/counts/uni_mapped_0err.bam',
                        'n_mapped_0err')
    yielded = []
    counts.array.each { |x| yielded << x }
    expect(yielded).to eq [['n_total', 625, 100],
                           ['n_mapped', 313, 50.080000000000005],
                           ['n_mapped_0err', 166, 26.56]]
  end

  it '.parse_nreads (expected failure)' do
    expect { counts.parse_nreads('./foobar', 'n_total') }
      .to raise_error(RuntimeError)
  end
end
