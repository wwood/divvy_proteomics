# divvy_proteomics

Takes a DTASelect CSV file, and parses the result so non-unique peptides get accounted for.

## Install
Get ruby somehow, if you don't already have it. Then, install this gem:
```
$ gem install divvy_proteomics
```
Or if there is permissions problems e.g. on OSX,
```
$ sudo gem install divvy_proteomics
```

## Usage
To test it work and to get a full listing of help
```sh
$ divvy_spectra -h
```

To run on a PepXML file e.g.
```
$ divvy_spectra --pep-xml my.pep.xml
```
or a DTASelect file
```
$ divvy_spectra DTASelect_file
```
Output is a table, with a row for each protein with a few columns, including number of unique spectra and the 
estimated number of spectral counts after sorting out the non-uniqueness.

Full usage information:
```

    Usage: divvy_spectra [options] <input_file>

    Takes a tab separated file containing a (possibly modified) output from a DTAselect run (or a pepXML file and add the flag --pep-xml), and use some algorithm to divy up the spectra that match multiple peptides.

        --merge-proteins FILE_OF_IDENTIFIERS
                                     Provide a space/tab separated file where the identifiers on each row should be treated as one protein
        --whitelist FILE_OF_PROTEINS_TO_REPORT
                                     Only report proteins that are in this whitelist, after divvying with everything
        --contaminant-regexes REGEXES
                                     Comma-separated list of regular expressions to apply to protein names. If the protein name matches then all spectra assigned to that protein are considered contaminants. [default: ]

Optional arguments:

        --pep-xml                    Input file is pep XML, rather than a DTA select output file [default: false]

Verbosity:

    -q, --quiet                      Run quietly, set logging to ERROR level [default INFO]
        --logger filename            Log to file [default stderr]
        --trace options              Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG
```

## Contributing to divvy\_proteomics

* Check out the latest master to make sure the feature hasn't been implemented or the bug hasn't been fixed yet.
* Check out the issue tracker to make sure someone already hasn't requested it and/or contributed it.
* Fork the project.
* Start a feature/bugfix branch.
* Commit and push until you are happy with your contribution.
* Make sure to add tests for it. This is important so I don't break it in a future version unintentionally.
* Please try not to mess with the Rakefile, version, or history. If you want to have your own version, or is otherwise necessary, that is fine, but please isolate to its own commit so I can cherry-pick around it.

## Copyright

Copyright (c) 2013-2015 Ben J Woodcroft. See LICENSE.txt for
further details.

