# bio-jaspar

## Tools for JASPAR motif analysis

[![Build Status](https://secure.travis-ci.org/wassermanlab/bioruby-jaspar.png)](http://travis-ci.org/wassermanlab/bioruby-jaspar)

This gem provides methods for: 

1. Reading and writing sequence motifs in JASPAR format
2. Accessing a JASPAR5 formatted database
3. Comparing, searching, and analyzing motifs in sequences

<sup>*</sup> **Note:** The JASPAR motif analysis tools consist of several modules that are directly imported from the Bio.motifs package in BioPython. Namely, those modules/submodules are: Bio.motifs, Bio.motifs.matrix, Bio.motifs.thresholds, Bio.motifs.jaspar. The functionality of this gem will be identical to the aforementioned modules/submodules.


## Installation

```sh
gem install bio-jaspar
```

## Usage

### Loading the gem

```ruby
require 'bio-jaspar'
```

### Loading a motif/motifs from a JASPAR database

A connection to the JASPAR database is made by creating a JASPAR5 instance.

```ruby
# Substitute the database credentials!
db = Bio::Jaspar::JASPAR5.new(
	:host => <db_host.org>,
	:name => <db_name>,
	:user => <db_user>,
	:password => <db_password>
)
```

Now, a motif can be retrieved by the matrix_id

```ruby
m = db.fetch_motif_by_id("MA0049")
puts m.to_s
```

Or multiple motifs can be retrieved by various criteria

```ruby
motifs = db.fetch_motifs(
	:collection => "CORE",
	:tax_group => ["fungi", "vertebrate"],
	:tf_class => "Helix-Turn-Helix",
	:min_ic => 2
)
motifs.each { |m| # do something with a motif }
```

### Motif analysis

Many methods are available for motif analysis. Here are some examples:

```ruby
m = db.fetch_motif_by_id("MA0049")

# Consensus sequence
m.consensus 					# BioRuby Sequence object
puts m.consensus

# Anticonsensus sequence
m.anticonsensus				# BioRuby Sequence object
puts m.anticonsensus

# Reverse complement motif
m.reverse_complement	# Bio::Motif::Motifs object

# Pseudocounts
m.pseudocounts

# Background
m.background

# Position weight matrix
m.pwm

# Position specific scoring matrix
m.pssm
```

Matrix methods are also available. Here are some examples:

```ruby
m = db.fetch_motif_by_id("MA0049")

# Maximum possible score for the given motif
m.pssm.max

# Minimum possible score for the given motif
m.pssm.min

# Expected value of the motif score
m.pssm.mean

# Standard deviation of the given motif score
m.pssm.std

# Find hits with the PWM score above given threshold
m.pssm.search(Bio::Sequence.auto("ACCTGCCTAAAAAA"), threshold = 0.5)
```

### Read/write Jaspar file

Already downloaded pfm, jaspar, sites files can be loaded/written using the Jaspar module

```ruby
# Read a pfm file
f = File.open("test.pfm", "r")
Bio::Jaspar.read(f, "pfm")
f.close

# Write motifs into a jaspar file
motifs = db.fetch_motifs(
	:collection => "CORE",
	:tax_group => ["fungi", "vertebrate"],
	:tf_class => "Helix-Turn-Helix",
	:min_ic => 2
)
File.open("test.jaspar", "w") do |f|
	Bio::Jaspar.write(f, "jaspar")
end
```

Please refer to the rdoc for full information on all available methods & classes.

## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

  http://github.com/wassermanlab/bioruby-jaspar

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

If you use this software, please cite one of
  
* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

## Biogems.info

This Biogem is published at (http://biogems.info/index.html#bio-jaspar)

## Copyright

Copyright (c) 2015 Jessica Lee. See LICENSE.txt for further details.

