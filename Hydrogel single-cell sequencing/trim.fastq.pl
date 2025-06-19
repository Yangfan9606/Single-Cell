#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-19
# version:   1.0
# license:   MIT
# brief:     Trim FASTQ of hydrogel-seq
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : 第二步，根据Barcode UMI check结果，提取可用reads，并trim (之后用fix.R1R2.clipper.pl做双端reads过滤,保留pair)

Usage: .pl [*BC_check.gz] [1.fastq.gz] [2.fastq.gz] -o [out.name].fastq.clipper.gz
        -o                      out name
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "help!");
##########
die $usage if ( @ARGV!=3 || defined($opts{"help"}));

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
my $infile2=shift;
open(INfile2,"zcat $infile2|") or die $!;
my $outname = $opts{o};
die "Input out name use -o\n" if $opts{o} eq "";

open(OUT1, "|gzip >$outname.1.fastq.clipper.gz") or die $!;
open(OUT2, "|gzip >$outname.2.fastq.clipper.gz") or die $!;

#############
my %Lake;
while(<INfile>){
        chomp;
        my @temp = split/\t/;
        my $n = $temp[0];
        my $b3 = $temp[2];
        my $b3_2 = $temp[3];
        my $b2 = $temp[4];
        my $b1 = $temp[5];
        my $len = length $temp[1];
        $len -= 5;
        if ($b3 ne "NA" && $b3 !~/not/ && $b3_2 ne "NA" && $b2 ne "NA" && $b2 !~/not/ && $b1 ne "NA" && $b1 !~/not/){
                $Lake{$n}{len} = 92 + $len;
        }
}
close INfile;
print "check read done...\n";
my %R1;
my $R1_tn5 = "CTGTCTCTTATACACA";
my %R1_b;
while(<INfile1>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile1>;
        my $L3 = <INfile1>;
        my $L4 = <INfile1>;
        chomp $L1;
        chomp $L2;
        chomp $L3;
        chomp $L4;
        my @temp = split/ /,$L1;
        my $n = $temp[0];
        $n =~ s/@//;
        next if $Lake{$n}{len} eq "";
        my $l = $Lake{$n}{len};
        my $bseq = substr($L2,0,$l);
        my @S = split/$bseq|$R1_tn5/,$L2;
#       print "@S, $S[0], $S[1], $S[2]\n";
        my $ok = $S[1];
        next if $ok eq "";
        $R1{$n} = 1;
        my $lok = length $ok;
        my $q = substr($L4,$l,$lok);
        my $r1_bc = substr($bseq,-15);
        my $vr1_bc = reverse $r1_bc;
        $vr1_bc =~ tr/ATGC/TACG/;
        $R1_b{$n} = $vr1_bc;
        print OUT1 "$L1\n$ok\n$L3\n$q\n";
}
close INfile1;
my %R2;
while(<INfile2>){
        next unless $_ =~ /^@/;
        my $L1 = $_;
        my $L2 = <INfile2>;
        my $L3 = <INfile2>;
        my $L4 = <INfile2>;
        chomp $L1;
        chomp $L2;
        chomp $L3;
        chomp $L4;
#       next if $L2 eq "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
        my @temp = split/ /,$L1;
        my $n = $temp[0];
        $n =~ s/@//;
        next if $R1{$n} eq "";
        my @S = split/$R1_b{$n}/,$L2;
        my $seq = $S[0];
        next if $seq eq "";
        if (@S > 1){
                my $l = length $seq;
                my $q = substr($L4,0,$l);
                print OUT2 "$L1\n$seq\n$L3\n$q\n";
        }else{
                print OUT2 "$L1\n$L2\n$L3\n$L4\n";
        }
}
close INfile2;

#############
close OUT;

Ptime("End");
print "##########End############\n";
