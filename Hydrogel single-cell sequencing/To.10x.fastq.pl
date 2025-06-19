#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-19
# version:   1.0
# license:   MIT
# brief:     Transfer FASTQ to 10x input format
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

Usage: .pl [10x-barcode.txt] [BAG600_2.2.fastq.clipper.gz] [BAG600_2.R1.BC_check.gz] [BAG600_2.1.fix.fastq.gz] -o [out.name]
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "help!");
##########
die $usage if ( @ARGV!=4 || defined($opts{"help"}));

###################################################
#                  START                          #
###################################################
my $optko;
foreach my $opt(keys %opts){
        $optko .= " -$opt $opts{$opt}";
}
print "##########Start############ perl $0 @ARGV ($optko)\n";
Ptime("Start");
my $barcodefile=shift;
open(BAR,"zcat $barcodefile|") or die $!;
my $infile=shift;
open(INfile,"zcat $infile|") or die $!;
my $infile1=shift;
open(INfile1,"zcat $infile1|") or die $!;
my $infile2=shift;
open(INfile2,"zcat $infile2|") or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT, "|gzip >$outname.For10x.1.fastq.gz") or die $!;
open(OUT2, "|gzip >$outname.For10x.2.fastq.gz") or die $!;

#############
my @BC10;
while(<BAR>){
        chomp;
        push @BC10,$_;
}
close BAR;
my %Lake;

while(<INfile>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile>;
        my $L3 = <INfile>;
        my $L4 = <INfile>;
        chomp $L1;
        chomp $L2;
        chomp $L3;
        chomp $L4;
        my @temp = split/ /,$L1;
        my $name = $temp[0];
        $name =~ s/^@//;
        $Lake{$name} = 1;
        my $R2seq = reverse $L2;
        $R2seq =~ tr/ATGC/TACG/;
        my $Q = reverse $L4;
        print OUT2 "$L1\n$R2seq\n$L3\n$Q\n";
}
close INfile;
my (%BC,%umi,%FBC);
while(<INfile1>){
        chomp;
        my @temp = split/\t/;
        my $name = $temp[0];
        next unless $Lake{$name} == 1;
        my $barcode = $temp[2].$temp[4].$temp[5];
        $barcode = substr($barcode,1,16);
        $BC{$name}=$barcode;
        $FBC{$barcode} = 1;
        $umi{$name}="AAAA$temp[3]$temp[-1]";
}
close INfile1;
my %TB;
my $l=0;
foreach my $b(keys %FBC){
        $TB{$b} = $BC10[$l];
        $l++;
}
while(<INfile2>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile2>;
        my $L3 = <INfile2>;
        my $L4 = "FFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        chomp $L1;
        chomp $L2;
        chomp $L3;
        my @temp = split/ /,$L1;
        my $name = $temp[0];
        $name =~ s/^@//;
        next unless $Lake{$name} == 1;
        my $seq = $TB{$BC{$name}}.$umi{$name};
        print OUT "$L1\n$seq\n$L3\n$L4\n";
}
close INfile;
#############
close OUT;

Ptime("End");
print "##########End############\n";
