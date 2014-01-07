require 'bio-logger'
Bio::Log::LoggerPlus.new('divvy_proteomics')
module Bio
  module DivvyProteomics
    module Logging
      def log
        Bio::Log::LoggerPlus['divvy_proteomics']
      end
    end
  end
end

require 'divvyable_protein'
require 'dta_select_output'
require 'pep_xml'

