require 'fileutils'
require 'spec_helper'
require 'ribopip/gene_id_2_name'

describe Ribopip::GeneID2Name do
  let(:id2name) {
    Ribopip::GeneID2Name.new('./spec/testdata/ref/mini.GRCm38.80.gtf')
  }

  it '.build_index' do
    expect(
      FileUtils.compare_file(
        './spec/testdata/ref/mini.GRCm38.80.gtf.name_index',
        './spec/testdata/ref/expected.name_index'
      )
    ).to be_truthy
  end

  it '.read_index' do
    foo2name = Ribopip::GeneID2Name.new('./spec/testdata/ref/expected')
    expect(foo2name.id2name).to eq id2name.id2name
  end

  it '.translate' do
    id2name.translate(
      1,
      './spec/testdata/gene_id_2_name/d1_vs_d2.deseq1.all.ft.csv',
      './spec/testdata/gene_id_2_name/d1_vs_d2.deseq1.names.csv'
    )

    expect(
      FileUtils.compare_file(
        './spec/testdata/gene_id_2_name/d1_vs_d2.deseq1.names.csv',
        './spec/testdata/gene_id_2_name/expected.names.csv'
      )
    ).to be_truthy
  end
end
