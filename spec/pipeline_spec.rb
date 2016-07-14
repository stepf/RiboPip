require 'spec_helper'
require 'ribopip/pipeline'

describe Ribopip::Pipeline::PipelineStep do
  let(:pipeline) { Ribopip::Pipeline::PipelineStep.new('foo', false) }
  it '.skip_step?: true; file exists' do
    expect(
      pipeline.skip_step?('./spec/testdata/SRR315601_625.fastq', 'foo')
    ).to be true
  end

  it '.skip_step?: false; @force_overwrite == true' do
    pipe = Ribopip::Pipeline::PipelineStep.new('foo', true)
    expect(
      pipe.skip_step?('./spec/testdata/SRR315601_625.fastq', 'foo')
    ).to be false
  end

  it '.skip_step?: false; file does not exists' do
    expect(
      pipeline.skip_step?('notexistant', 'foo')
    ).to be false
  end
end
