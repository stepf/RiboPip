require 'fileutils'
require 'spec_helper'
require 'ribopip/expression_analysis'
require 'ribopip/deseq'
require 'ribopip/filenames'

names = Ribopip::Filenames::Deseq.new(
  'd1', 'd2', './spec/testdata/expression_analysis'
)
names_rdiff = Ribopip::Filenames::Ribodiff.new(
  'd1', 'd2', './spec/testdata/expression_analysis'
)

describe Ribopip::ExpressionAnalysis::Deseq::Deseq do
  let(:deseq) do
    Ribopip::ExpressionAnalysis::Deseq::Deseq.new(
      names.dup,
      [['./spec/testdata/expression_analysis/dummy1_rep1.txt',
        './spec/testdata/expression_analysis/dummy1_rep2.txt'],
       ['./spec/testdata/expression_analysis/dummy2_rep1.txt',
        './spec/testdata/expression_analysis/dummy2_rep2.txt']],
      %w(d1 d2),
      './spec/testdata/expression_analysis'
    )
  end

  it '.read_featurecounts' do
    deseq.read_featurecounts(
      './spec/testdata/expression_analysis/dummy1_rep1.txt', 0
    )
    expect(deseq.features).to eq(
      ENSMUSG00000051951: [59],
      ENSMUSG00000102851: [1],
      ENSMUSG00000033845: [18],
      ENSMUSG00000025903: [1],
      ENSMUSG00000033813: [14],
      ENSMUSG00000002459: [7],
      ENSMUSG00000033793: [135],
      ENSMUSG00000025905: [1],
      ENSMUSG00000033774: [1],
      ENSMUSG00000025907: [20],
      noname: [100]
    )
  end

  it '.merge_featurecounts' do
    deseq.merge_featurecounts
    expect(deseq.features). to eq(
      ENSMUSG00000002459: [7, 6, 2, 2],
      ENSMUSG00000025903: [1, 2, 1],
      ENSMUSG00000025905: [1],
      ENSMUSG00000025907: [20, nil, 14, 9],
      ENSMUSG00000033774: [1],
      ENSMUSG00000033793: [135, 230, 53, 42],
      ENSMUSG00000033813: [14, 26, 3, 1],
      ENSMUSG00000033845: [18, 39, 4, 8],
      ENSMUSG00000051951: [59, 129, 21, 21],
      noname: [100, 75, 50, 25],
      ENSMUSG00000102851: [1],
      ENSMUSG00000025902: [nil, 1],
      ENSMUSG00000103201: [nil, 1],
      ENSMUSG00000103922: [nil, 2, nil, 1],
      ENSMUSG00000104352: [nil, 1],
      ENSMUSG00000025909: [nil, nil, 2],
      ENSMUSG00000051285: [nil, nil, 10],
      ENSMUSG00000061024: [nil, nil, 5],
      ENSMUSG00000033740: [nil, nil, nil, 2],
      ENSMUSG00000102653: [nil, nil, nil, 1],
      ENSMUSG00000103845: [nil, nil, nil, 1]
    )
  end

  it '.write_merged_featurecounts' do
    deseq.merge_featurecounts
    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/d1_vs_d2.counts',
        './spec/testdata/expression_analysis/expected_d1_vs_d2.counts'
      )
    ).to be_truthy
  end
end

describe Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults do
  let(:add_counts) do
    Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
      './spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.tsv'
    )
  end
  let(:add_norm_counts) do
    Ribopip::ExpressionAnalysis::Deseq::AddCountsToResults.new(
      './spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.tsv',
      '_norm',
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.log'
    )
  end

  it '.parse_counts' do
    add_counts.parse_counts(
      './spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.tsv'
    )
    expect(add_counts.features[:noname]).to eq([100, 75, 50, 25])
  end

  it '.parse_rep_names' do
    names = add_counts.parse_rep_names(
      './spec/testdata/expression_analysis/d1_vs_d2.deseq.counts.tsv'
    )
    expect(names).to eq(%w(d1_rep1 d1_rep2 d2_rep1 d2_rep2))
  end

  it '.parse_norm_factors' do
    norm_factors = add_norm_counts.parse_norm_factors(
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.log'
    )
    expect(norm_factors).to eq([2, 2, 2, 2])
  end

  it '.compute (raw counts)' do
    add_counts.compute(
      1,
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.csv',
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.csv'
    )
    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.csv',
        './spec/testdata/expression_analysis/expected.all.ft.csv'
      )
    ).to be_truthy
  end

  it '.compute (with normalized counts & suffix)' do
    add_norm_counts.compute(
      1,
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.csv',
      './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.norm.csv'
    )
    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/d1_vs_d2.deseq1.all.counts.norm.csv',
        './spec/testdata/expression_analysis/expected.all.ft.norm.csv'
      )
    ).to be_truthy
  end
end

describe Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel do
  it '.convert' do
    Ribopip::ExpressionAnalysis::Deseq::Deseq2Excel.convert(
      './spec/testdata/expression_analysis/names.csv',
      './spec/testdata/expression_analysis/names.xml'
    )

    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/names.xml',
        './spec/testdata/expression_analysis/expected_names.xml'
      )
    ).to be_truthy
  end
end

describe Ribopip::ExpressionAnalysis::Ribodiff do
  let(:rdiff) do
    Ribopip::ExpressionAnalysis::Ribodiff.new(
      names_rdiff.dup,
      [['./spec/testdata/expression_analysis/dummy1_rep1.txt',
        './spec/testdata/expression_analysis/dummy1_rep2.txt',
        './spec/testdata/expression_analysis/dummy1_rep3.txt'],
       ['./spec/testdata/expression_analysis/dummy2_rep1.txt',
        './spec/testdata/expression_analysis/dummy2_rep2.txt'],
       ['./spec/testdata/expression_analysis/dummy1_rep1.txt',
        './spec/testdata/expression_analysis/dummy1_rep2.txt',
        './spec/testdata/expression_analysis/dummy1_rep3.txt'],
       ['./spec/testdata/expression_analysis/dummy2_rep1.txt',
        './spec/testdata/expression_analysis/dummy2_rep2.txt']],
      %w(fp_1 fp_2 mrna_1 mrna_2),
      './spec/testdata/expression_analysis'
    )
  end

  it '.compute_exp_outline' do
    rdiff.compute_exp_outline

    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/d1_vs_d2.rdiff.outline.csv',
        './spec/testdata/expression_analysis/expected.rdiff.outline.csv'
      )
    ).to be_truthy
  end

  it '.filter_results' do
    rdiff.filter_results

    expect(
      FileUtils.compare_file(
        './spec/testdata/expression_analysis/d1_vs_d2.rdiff.results.padj.tsv',
        './spec/testdata/expression_analysis/expected.rdiff.results.padj.tsv'
      )
    ).to be_truthy
  end
end
