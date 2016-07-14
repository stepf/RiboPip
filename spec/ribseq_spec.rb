require 'ribopip'

describe Ribopip do
  subject { Object.new.extend(Ribopip) }
  it '.run_cmd' do
    expect(subject.run_cmd('ls')).to eq true
  end

  it '.run_cmd (expected failure)' do
    expect { subject.run_cmd('ls non_existent') }.to raise_error(SystemCallError)
  end
end
