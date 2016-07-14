require 'spec_helper'
require 'ribopip/binary_checker'
describe Ribopip::BinaryChecker do
  context 'empty BinaryChecker object' do
    let(:bincheck) { Ribopip::BinaryChecker.new([]) }
    it '.get_binary_path' do
      path_expected = `which ruby`
      path_found = bincheck.get_binary_path('ruby')
      expect(path_found).to eq path_expected
    end

    it '.get_binary_version' do
      version_expected = RUBY_VERSION
      version_found = bincheck.get_binary_version('ruby', '-v', 'ruby')
      expect(version_found).to eq version_expected
    end

    it '.compare_versions: inequal' do
      warning = bincheck.compare_versions('2.0.0', '2.0.1')
      expect(warning).to be_a String
    end

    it '.compare_versions: equal' do
      warning = bincheck.compare_versions('2.0.0', '2.0.0')
      expect(warning).to eq nil
    end

    it '.check_binary: expected failure' do
      expect {
        bincheck.check_binary('foobarxyz', '-flag', 'version', 'regex')
      }.to raise_error(SystemExit)
    end

    it '.check_binary: expected sucess' do
      expect {
        bincheck.check_binary('ruby', '-v', "#{RUBY_VERSION}", 'ruby')
      }.to_not raise_error
    end
  end

  it '.perform' do
    expect {
      bin_array = [
        ['ruby', '-v', "#{RUBY_VERSION}", 'ruby'],
        ['ruby', '-v', "#{RUBY_VERSION}", 'ruby']
      ]
      bin_check = Ribopip::BinaryChecker.new(bin_array)
      bin_check.perform_checks
    }.to_not raise_error
  end
end
