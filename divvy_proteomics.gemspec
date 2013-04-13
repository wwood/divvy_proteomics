# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "divvy_proteomics"
  s.version = "0.0.1"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Ben J Woodcroft"]
  s.date = "2013-04-13"
  s.description = "divvy up spectra from DTASelect files in a somewhat parsimonious way"
  s.email = "donttrustben@gmail.com"
  s.executables = ["divvy_spectra"]
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.md"
  ]
  s.files = [
    ".document",
    ".rspec",
    "Gemfile",
    "LICENSE.txt",
    "README.md",
    "Rakefile",
    "VERSION",
    "bin/divvy_spectra",
    "divvy_proteomics.gemspec",
    "lib/divvy_proteomics.rb",
    "spec/data/merge_definition.csv",
    "spec/data/multiply_mapped_spectra.csv",
    "spec/data/single_protein.csv",
    "spec/data/single_protein_with_aliases.csv",
    "spec/data/three_proteins.csv",
    "spec/data/three_proteins_meant_for_merge.csv",
    "spec/data/three_proteins_with_contaminant.csv",
    "spec/divvy_proteomics_spec.rb",
    "spec/spec_helper.rb"
  ]
  s.homepage = "http://github.com/wwood/divvy_proteomics"
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = "1.8.24"
  s.summary = "divvy up spectra from DTASelect files in a parsimonious way"

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio-logger>, [">= 0"])
      s.add_development_dependency(%q<systemu>, [">= 0"])
      s.add_development_dependency(%q<rspec>, [">= 2.8.0"])
      s.add_development_dependency(%q<rdoc>, [">= 3.12"])
      s.add_development_dependency(%q<bundler>, [">= 1.0.0"])
      s.add_development_dependency(%q<jeweler>, [">= 1.8.4"])
    else
      s.add_dependency(%q<bio-logger>, [">= 0"])
      s.add_dependency(%q<systemu>, [">= 0"])
      s.add_dependency(%q<rspec>, [">= 2.8.0"])
      s.add_dependency(%q<rdoc>, [">= 3.12"])
      s.add_dependency(%q<bundler>, [">= 1.0.0"])
      s.add_dependency(%q<jeweler>, [">= 1.8.4"])
    end
  else
    s.add_dependency(%q<bio-logger>, [">= 0"])
    s.add_dependency(%q<systemu>, [">= 0"])
    s.add_dependency(%q<rspec>, [">= 2.8.0"])
    s.add_dependency(%q<rdoc>, [">= 3.12"])
    s.add_dependency(%q<bundler>, [">= 1.0.0"])
    s.add_dependency(%q<jeweler>, [">= 1.8.4"])
  end
end
