#!/usr/bin/perl
use strict;
use warnings;

#find cg-dinucleotides in gene sequences

my($in_file, $out_file) = @ARGV;

open GENE, "<$in_file";
chomp(my @lines = <GENE>);
open OUT, ">$out_file";


my $string = join('',@lines);
my @seqs = split(/>/,$string);


my %cg_count;
my %tg_count;
my %ca_count;
my %gc_count;
my $gene;
my $sequence;
my %c;
my %g;
my %t;
my %a;
my %length;  

foreach my $seq(@seqs){
   if($seq =~ /^(\S+)\s+\S+\d([ACTGN]+)$/i){
# $1: gene name; $2 sequence
      $gene = $1;
      $sequence = $2;
      @{$cg_count{$gene}} = $sequence =~ /CG/g;
      @{$tg_count{$gene}} = $sequence =~ /TG/g;
      @{$ca_count{$gene}} = $sequence =~ /CA/g;
      @{$gc_count{$gene}} = $sequence =~ /GC/g;
      @{$c{$gene}} = $sequence =~ /C/g;
      @{$g{$gene}} = $sequence =~ /G/g;
      @{$a{$gene}} = $sequence =~ /A/g;
      @{$t{$gene}} = $sequence =~ /T/g;
      $length{$gene} = length($sequence);
   }
   elsif($seq =~ /^(\S+\d)([ACTGN]+)$/i){
# $1: gene name; $2 sequence
      $gene = $1;
      $sequence = $2;
      @{$cg_count{$gene}} = $sequence =~ /CG/g;
      @{$tg_count{$gene}} = $sequence =~ /TG/g;
      @{$ca_count{$gene}} = $sequence =~ /CA/g;
      @{$gc_count{$gene}} = $sequence =~ /GC/g;
      @{$c{$gene}} = $sequence =~ /C/g;
      @{$g{$gene}} = $sequence =~ /G/g;
      @{$a{$gene}} = $sequence =~ /A/g;
      @{$t{$gene}} = $sequence =~ /T/g;
      $length{$gene} = length($sequence);
   }
}

my $cg_count;
my $tg_count;
my $ca_count;
my $gc_count;
my $c_count;
my $g_count;
my $t_count;
my $a_count;
my $cg_oe;
my $tg_oe;
my $ca_oe;
my $gc_oe;


print OUT "Gene\tCG\tTG\tCA\tGC\tC\tG\tT\tA\tlength\tCpG_oe\tTpG_oe\tCpA_oe\tGpC_oe\n";

foreach my $key (sort keys %cg_count){
   $cg_count = scalar(@{$cg_count{$key}});
   $tg_count = scalar(@{$tg_count{$key}});
   $ca_count = scalar(@{$ca_count{$key}});
   $gc_count = scalar(@{$gc_count{$key}});
   $c_count = scalar(@{$c{$key}});
   $g_count = scalar(@{$g{$key}});
   $t_count = scalar(@{$t{$key}});
   $a_count = scalar(@{$a{$key}});

if(($c_count * $g_count)/$length{$key} == 0){
   $cg_oe = 0;
   $gc_oe = 0;
}
else{
   $cg_oe = $cg_count / (($c_count * $g_count)/$length{$key});
   $gc_oe = $gc_count / (($c_count * $g_count)/$length{$key});
}

if(($t_count * $g_count)/$length{$key} == 0){
   $tg_oe = 0;
}
else{
   $tg_oe = $tg_count / (($t_count * $g_count)/$length{$key});
}

if(($c_count * $a_count)/$length{$key} == 0){
   $ca_oe = 0;
}
else{
   $ca_oe = $ca_count / (($c_count * $a_count)/$length{$key});
}

   print OUT "$key\t$cg_count\t$tg_count\t$ca_count\t$gc_count\t$c_count\t$g_count\t$t_count\t$a_count\t$length{$key}\t$cg_oe\t$tg_oe\t$ca_oe\t$gc_oe\n";
}

