#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $filename = $ARGV[0];

open IN, $filename;

my %samples;

while(my $line = <IN>) {
  chomp($line);
  my @x =split(/\s/,$line);
  $samples{$x[0]}=$x[1];
}

my @k_sample = keys %samples;

my $header = "";
my $seq = "";
my %proteins;

foreach my $sample (@k_sample ) {
  my $file = $samples{$sample};
  open IN, $file;
  while (my $line = <IN>) {
    chomp($line);
    if ($line =~ /^>/) {
      $header = $line;
      $header =~ s/>//;
    }
    else {
      $seq = $line;
      $proteins{$header}{$sample}.=$seq;
    }
  }
  close IN;
}

my %final_prots;

my @l_sample = keys %samples;

foreach my $protein (keys %proteins) {
  my @check_diffs;
  foreach my $sample (@k_sample) {
    push (@check_diffs,$proteins{$protein}{$sample});
  }
  my @uniq_prots = uniq (@check_diffs);
  if (scalar(@uniq_prots)>1 && length($proteins{$protein}{$k_sample[0]})%3==0) {
    open OUT, ">protein_fastas/$protein".".aln";
    #open OUT, ">test";
    foreach my $sample (@k_sample) {
      $final_prots{$sample}.= $proteins{$protein}{$sample};
      print OUT ">$sample\n$proteins{$protein}{$sample}\n";
    }
    close OUT;
  }
}

#foreach my $key (keys %final_prots) {
#  print ">$key\n$final_prots{$key}\n";
#}
