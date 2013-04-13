# divvy_proteomics

Takes a DTASelect CSV file, and parses the result so non-unique peptides get accounted for.

## Install
Get ruby somehow, if you don't already have it. Then, install this gem:
```
gem install divvy_spectra
```

## Usage
```
$ divvy_spectra <DTASelectFile>
```
Output is a table, with a row for each protein with a few columns, including number of unique spectra and the 
estimated number of spectral counts after sorting out the non-uniqueness.

Full usage information:
```
$ divvy_spectra -h

    Usage: divvy_spectra [options] <DTASelect_file>

    Takes a tab separated file containing a (possibly modified) output from a DTAselect run, and use some algorithm to divy up the spectra that match multiple peptides.

        --merge-proteins FILE_OF_IDENTIFIERS
                                     Provide a space/tab separated file where the identifiers on each row should be treated as one protein
        --whitelist FILE_OF_PROTEINS_TO_REPORT
                                     Only report proteins that are in this whitelist, after divvying with everything

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

Copyright (c) 2013 Ben J Woodcroft. See LICENSE.txt for
further details.

