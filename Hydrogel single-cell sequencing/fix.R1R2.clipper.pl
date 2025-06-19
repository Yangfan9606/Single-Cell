#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-19
# version:   1.0
# license:   MIT
# brief:     Fix reads for cellranger input
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : Processing the file

Usage: .pl IN_file -o [out.name]
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "help!");
##########
die $usage if ( @ARGV!=2 || defined($opts{"help"}));

###################################################
#                  START                          #
###################################################
my $optko;
foreach my $opt(keys %opts){
        $optko .= " -$opt $opts{$opt}";
}
print "##########Start############ perl $0 @ARGV ($optko)\n";
Ptime("Start");
my $infile=shift;
open(INfile,"zcat $infile|") or die $!;
my $infile1=shift;
open(INfile1,"zcat $infile1|") or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT1, "|gzip >$outname.1.fix.fastq.gz") or die $!;
open(OUT2, "|gzip >$outname.2.fix.fastq.gz") or die $!;

#############
my %Lake;

while(<INfile>){
        next unless $_ =~ /^@/;
        chomp;
        my @temp = split/ /;
        my $L2 = <INfile>;
        chomp $L2;
        $Lake{$temp[0]} += 1 if $L2 ne "";
}
close INfile;
while(<INfile1>){
        next unless $_ =~ /^@/;
        chomp;
        my @temp = split/ /;
        my $L2 = <INfile1>;
        chomp $L2;
        $Lake{$temp[0]} += 1 if $L2 ne "";
}
close INfile1;
open(INfile,"zcat $infile|") or die $!;
while(<INfile>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile>;
        my $L3 = <INfile>;
        my $L4 = <INfile>;
        chomp $L1;
        my @temp = split/ /,$L1;
        next unless $Lake{$temp[0]} == 2;
        print OUT1 "$L1\n$L2$L3$L4";
}
close INfile;
open(INfile1,"zcat $infile1|") or die $!;
while(<INfile1>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile1>;
        my $L3 = <INfile1>;
        my $L4 = <INfile1>;
        chomp $L1;
        my @temp = split/ /,$L1;
        next unless $Lake{$temp[0]} == 2;
        print OUT2 "$L1\n$L2$L3$L4";
}
close INfile1;
#############
close OUT;

Ptime("End");
print "##########End############\n";
