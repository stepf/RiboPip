# version constant & helper functions
module Ribopip
  # version constant for https://github.com/mojombo/rakegem
  VERSION = '0.1.0'

  # Output message to $stderr, prefixed with the program name and timestamp
  #
  # *args - any number of strings
  def print_e(*args)
    timestamp = Time.now.strftime('%Y-%m-%d %H:%M:%S')
    args.each do |msg|
      STDERR.puts("[Ribopip #{timestamp}] #{msg}")
    end
  end

  # Runs system call verbosely, fails if exit code != 0
  #
  # cmd - cmd string to be run
  #
  # Returns boolean
  def run_cmd(cmd)
    print_e("Run: #{cmd}")
    exit_code = system(cmd)
    fail SystemCallError, cmd unless exit_code
    exit_code
  end
end
