require 'spec_helper'
require 'ribopip/filenames'

describe Ribopip::Filenames::Alignment do
  let(:names) { Ribopip::Filenames::Alignment.new('./dir/SRR315601_625.fastq') }
  it '.new' do
    expect(names.base).to eq './dir/SRR315601_625'
  end

  it '.get' do
    expect(names.get('counts')).to eq './dir/SRR315601_625.counts'
  end

  it '.set' do
    names.set('counts', '.countsfoo')
    expect(names.get('counts')).to eq './dir/SRR315601_625.countsfoo'
  end

  it '.set_bucket' do
    names.set_bucket(10, 20)
    expect(names.get('counts')).to eq './dir/10.20.SRR315601_625.counts'
  end

  it '.unset_bucket' do
    names.set_bucket(10, 20)
    expect(names.get('counts')).to eq './dir/10.20.SRR315601_625.counts'
    names.unset_bucket
    expect(names.get('counts')).to eq './dir/SRR315601_625.counts'
  end

  it '.get_info' do
   expect(names.get_info('counts')).to eq '#Reads after each alignment step, '\
                                          'inluding buckets'
  end
end
