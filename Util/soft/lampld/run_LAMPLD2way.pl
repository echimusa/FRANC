#!/usr/local/bin/perl -w
use strict;

if(@ARGV != 5){
    print "USAGE::perl run_LAMPLD.pl <posfile> <EURHaps> <YRIHaps> <genfile> <outfile>\n";
    die;
}

my $numStatesHMM = 50;
my $winSize = 50;
my $posfile = shift @ARGV;
my $EURfile = shift @ARGV;
my $YRIfile = shift @ARGV;
my $genfile = shift @ARGV;
my $outfile = shift @ARGV;

my $cmd ="./bin/unolanc2way $winSize $numStatesHMM $posfile $EURfile $YRIfile $genfile $outfile";
print "Running LAMPLD::$cmd\n";

`$cmd`;


