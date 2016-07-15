require 'fileutils'
require 'spec_helper'
require 'ribopip/metrics'

describe Ribopip::Metrics::GtSeqstat do
  it '.compute' do
    Ribopip::Metrics::GtSeqstat.compute(
      './spec/testdata/SRR315601_625.fastq',
      './spec/testdata/statistics/SRR315601_625.fastq.gt_seqstat'
    )
    # only the first 6 lines are relevant
    `head -n 6 ./spec/testdata/statistics/SRR315601_625.fastq.gt_seqstat > ./spec/testdata/statistics/SRR315601_625.fastq.gt_seqstat`
    expect(
      FileUtils.compare_file(
        './spec/testdata/statistics/SRR315601_625.fastq.gt_seqstat',
        './spec/testdata/statistics/expected.gt_seqstat'
      )
    ).to be_truthy
  end
end
