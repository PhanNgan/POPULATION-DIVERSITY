#!/usr/bin/env perl

use strict;
use warnings;


my $reference = $ARGV[0];
my $vcf = $ARGV[1];

my $header = "";
my $seq = "";
my %ref;

open IN, $reference;

while (my $line = <IN>) {
  chomp($line);
  if ($line =~ /^>/) {
    $header = $line;
    $header =~ s/>//;
  }
  else {
    $seq = $line;
    $ref{$header}.=$seq;
  }
}

close IN;

open IN,$vcf;
my $k = 0;
while (my $line=<IN>) {
  chomp($line);
  unless ($line =~ /^#/) {
    my @p = split(/\t/,$line);
    my $contig = $p[0];
    my $loc = $p[1]-1;
    my $ref_base = $p[3];
    my $alt_base = $p[4];
    my $genotype = $p[9];
    if ($genotype =~ /^1\/1/) {
      my @new_ref = split("",$ref{$contig});
      my @alts = split("",$alt_base);
      my $i = 0;
      while ($i<@alts) {
        $new_ref[$loc+$i]=$alts[$i];
        $i++;
      }
      my $new_reference = join("",@new_ref);
      $ref{$contig} = $new_reference;
    }
    $k++;
    if ($k%1000 == 0) {
      print STDERR "$k finished\n";
    }
  }
}

foreach my $key (keys %ref) {
  print ">$key\n$ref{$key}\n";
}
