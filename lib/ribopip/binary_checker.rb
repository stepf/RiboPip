module Ribopip
  # Check if all required binaries are present on the host system.
  # Exit if not. Warn if version is different. Input format: Array of 4-tupel:
  # [binary, flag, version, regex],
  # e.g. ['samtools', '--version', '1.1', 'regex']
  class BinaryChecker
    include Ribopip

    # Initializes array
    #
    # check_array - array containing [binary, flag, version, regex]
    def initialize(check_array)
      @check_array = check_array
    end

    # Iterates over array and runs check_binary on every four-tupel.
    def perform_checks
      @check_array.each do |binary, flag, version, regex|
        check_binary(binary, flag, version, regex)
      end
    end

    # Generic wrapper function to get and print binary path and version number.
    # Aborts if binary is not found.
    #
    # binary  - binary name
    # flag    - version flag
    # version - expected version
    # regex   - string used for extraction
    #
    # Example
    #   check_binary('samtools', '--version', '1.1', 'samtools')
    def check_binary(binary, flag, version, regex)
      path = get_binary_path(binary)
      if path == ''
        abort "#{binary} not found. Please make sure it is in " \
              'properly installed and your PATH correctly set' \
              'and start pipeline again.'
      end

      version_found = get_binary_version(binary, flag, regex)
      comparison = compare_versions(version, version_found)

      print_e("Using #{binary} #{version_found}: #{path}")
      print_e("#{comparison}") unless comparison.nil?
    end

    # Gets path to binary using system command.
    #
    # binary - name of binary
    #
    # Returns string.
    def get_binary_path(binary)
      `which #{binary}`
    end

    # Gets binary version by calling it with flag and extracting regex
    #
    # binary - name of binary
    # flag   - version flag
    # regex  - string used for extraction
    #
    # Example
    #   get_binary_version('samtools', '--version', 'samtools')
    #
    # Returns string or nil (if regex did not hit).
    def get_binary_version(binary, flag, regex)
      `#{binary} #{flag} 2>&1`[/#{regex}\s*[\d\.]+/][/[\d\.]+/] rescue 'NA'
    end

    # Prints warning if binary version differs from recommended one.
    #
    # version       - recommended version
    # version_found - found binary version
    #
    # Returns nothing.
    def compare_versions(version, version_found)
      # Split version numbers into arrays to make multiple dots comparable.
      "WARNING: #{version_found} is different from #{version} (recommended). " \
      'Results might deviate.' if version != version_found
    end
  end
end
