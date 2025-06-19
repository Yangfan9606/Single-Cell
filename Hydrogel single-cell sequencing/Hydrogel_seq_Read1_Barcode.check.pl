#!/usr/bin/perl
#------------------------------------------#
# author:    Yangfan Zhou
# email:     yangfanzhou9606@gmail.com
# date:      2025-06-19
# version:   1.0
# license:   MIT
# brief:     Check barcode in FASTQ
#------------------------------------------#
use strict;
use Getopt::Long;
use lib './';
use yangfan;

my %opts;
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program : 第一步，检查R1的barcode 分布;INPUT Reads1.fastq (name.fastq.gz)

Usage: .pl R1.FASTQ_file -o [out.name].R1.BC_check
        -o                      out name
        -b1                     barcode1 file (TSO 3') [barcode1_TSO.txt]
        -b2                     barcode2 file [barcode2.txt]
        -b3                     barcode3 file (5') [barcode3.txt]
        -help                   output help information

USAGE

GetOptions(\%opts, "l:s", "b:s", "o:s", "u:s", "i:s", "b1:s", "b2:s", "b3:s", "help!");
##########
die $usage if ( @ARGV!=1 || defined($opts{"help"}));

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
open(INfile, "zcat $infile|") or die $!;
my $outname = $opts{o};
my (%B1,%B2,%B3);
my $b1f = $opts{b1};
open(IB1,$b1f) or die $!;
my $b2f = $opts{b2};
open(IB2,$b2f) or die $!;
my $b3f = $opts{b3};
open(IB3,$b3f) or die $!;
die "Input out name use -o\n" if $opts{o} eq "";
open(OUT, "|gzip>$outname.R1.BC_check.gz") or die $!;
while(<IB1>){
        chomp;
        $B1{$_} = 1;
}
close IB1;
while(<IB2>){
        chomp;
        $B2{$_} = 1;
}
close IB2;
while(<IB3>){
        chomp;
        $B3{$_} = 1;
}
close IB3;
#############
my (%Lake,%name);
my ($del,$ins,$ok);
my $line = 0;
my $TU = "AGTGGAAAAGGAAGGTGGT";
my $TU_mis1 = fuzzy_pattern($TU,1);
my $len_TU = length $TU;
my $b2seq = "GACCGACTCGCATTACCCTAT";
my $b2seq_mis1 = fuzzy_pattern($b2seq,1);
my $len_b2s = length $b2seq;
my $b3seq = "TCTTGACACGAAGGATAG";
my $b3seq_mis1 = fuzzy_pattern($b3seq,1);
while(<INfile>){
        next unless $_ =~ /^@/;
        $line++;
        my $L1 = $_;
        my $L2 = <INfile>;
        my $L3 = <INfile>;
        my $L4 = <INfile>;
        chomp $L1;
        chomp $L2;
        chomp $L3;
        chomp $L4;
        my @n = split/ /,$L1;
        $n[0] =~ s/@//;
        my $name = $n[0];
        my $seq = $L2;
        my @S1 = split/$TU_mis1/,$seq;
        my ($r1,$r2,$r3,$r4,$r5,$r6,$r7)=("NA","NA","NA","NA","NA","NA","NA");
        my $bag_UMI;
########### 5' , b3, 2N, b2, b1,GGG,len
        my $r8;
        if (@S1 == 2){
#my $TU = "AGTGGAAAAGGAAGGTGGT";
#my $b2seq = "GACCGACTCGCATTACCCTAT";
#my $b3seq = "TCTTGACACGAAGGATAG";
                my @S2 = split/$b2seq_mis1/,$S1[0];
                if (@S2 == 2){
                        my @S3 = split/$b3seq_mis1/,$S2[0];
                        if (@S3 == 2){
                                $r1 = $S3[0];
                                my $lS31 = length $S3[1];
                                if ($lS31 != 8){
                                        $r2 = "not8";
                                        $r3 = "not8";
                                }else{
                                        my @base = split//,$S3[1];
                                        $r3 = $base[-2].$base[-1];
                                        pop @base;
                                        pop @base;
                                        my $rB3 = join("",@base);
                                        $r2 = $B3{$rB3} eq ""?"NA":$rB3;
                                }
                        }
                        my $lS21 = length $S2[1];
                        if ($lS21 != 6){
                                $r4 = "not6";
                        }else{
                                $r4 = $B2{$S2[1]} eq ""?"NA":$S2[1];
                        }
                }else{
                        $r4 = "NA";
                        my @S3 = split/$b3seq_mis1/,$S1[0];
                        if (@S3 == 2){
                                $r1 = $S3[0];
                                my $S3_1 = substr($S3[1],6);
                                my @base = split//,$S3_1;
                                $r3 = $base[-2].$base[-1];
                                pop @base;
                                pop @base;
                                my $rB3 = join("",@base);
                                $r2 = $B3{$rB3} eq ""?"NA":$rB3;
                        }else{
                                $r1="NA";
                                $r2="NA";
                                $r3="NA";
                        }
                }
#my $TU = "AGTGGAAAAGGAAGGTGGT";
#my $b2seq = "GACCGACTCGCATTACCCTAT";
#my $b3seq = "TCTTGACACGAAGGATAG";
                my @TSO_umi = split//,$S1[1];
                my $umi_1 = $TSO_umi[0].$TSO_umi[1].$TSO_umi[2].$TSO_umi[3].$TSO_umi[4].$TSO_umi[5];
                $bag_UMI = $TSO_umi[6].$TSO_umi[7].$TSO_umi[8].$TSO_umi[9].$TSO_umi[10].$TSO_umi[11];
                my $GGG = $TSO_umi[12].$TSO_umi[13].$TSO_umi[14];
                $r5 = $B1{$umi_1} eq ""?"NA":$umi_1;
                if ($GGG =~ /GGG$/ || $GGG =~ /[ATGC]GG$/ || $GGG =~ /G[ATGC]G$/ || $GGG =~ /GG[ATGC]$/){
                        $r6 = $GGG;
                }else{
                        $r6 = "NA";
                }
        }else{
                $r5 = "NA";
                $r6 = "NA";
#my $b2seq = "GACCGACTCGCATTACCCTAT";
#my $b3seq = "TCTTGACACGAAGGATAG";
                my @S2 = split/$b2seq/,$seq;
                if (@S2 == 2){
                        my @S3 = split/$b3seq_mis1/,$S2[0];
                        $r1 = $S3[0];
                        my $lS31 = length $S3[1];
                        if ($lS31 != 8){
                                $r2 = "not8";
                                $r3 = "not8";
                        }else{
                                my @base = split//,$S3[1];
                                $r3 = $base[-2].$base[-1];
                                pop @base;
                                pop @base;
                                my $rB3 = join("",@base);
                                $r2 = $B3{$rB3} eq ""?"NA":$rB3;
                        }
                }else{
                        my @S3 = split/$b3seq_mis1/,$seq;
                        if (@S3 == 2){
                                $r1 = $S3[0];
                                my $S3_1 = substr($S3[1],6);
                                my @base = split//,$S3_1;
                                $r3 = $base[-2].$base[-1];
                                pop @base;
                                pop @base;
                                my $rB3 = join("",@base);
                                $r2 = $B3{$rB3} eq ""?"NA":$rB3;
                        }else{
                                $r1="NA";
                                $r2="NA";
                                $r3="NA";
                        }
                }
        }
        if ($r1 ne "NA" && $r2 !~/not/ && $r2 ne "NA" && $r4 !~/not/ && $r4 ne "NA" && $r5 ne "NA" && $r6 ne "NA"){
                my $l1 = length $r1;
#               $r7 = $l1 + 8 + 6 + 6 + 3 + 19 + 21 + 18 + 6;
                $r8 = "OK";
        }else{
#               $r7 = "NA";
                $r8 = "notOK";
        }
        print OUT "$name\t$r1\t$r2\t$r3\t$r4\t$r5\t$r6\t$r8\t$bag_UMI\n";
        print "$line reads prcessed ... \n" if (($line % 1000000) == 0);
}
close INfile;
#############
close OUT;

Ptime("End");
print "##########End############\n";
