require 'bio-logger'

module Bio
  module DivvyProteomics
    module Logging
      def log
        Bio::Log::LoggerPlus['divvy_proteomics']
      end
    end
  end
end

require 'dta_select_output'
require 'pep_xml'
