#!/usr/bin/perl
use strict;
use warnings;

###filter the FORMAT/AD field of a vcf file that has been called using the following bcftools (v1.8) command:
#bcftools mpileup -b list_of_bam_or_cram_files.txt -R genomic_regions_to_call.bed -f reference_file.fa -C 0 -q 20 -Q 20 -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR -Ou | bcftools call -Oz --ploidy 1 -m -o outout.vcf.gz


####INPUT and OUTPUT files - specify on the command like with full paths
#$input - input vcf
#$output - output vcf
#$output_changed_lines - output file that lists all the changed genotypes
#$multi_changed_lines - output file that lists all multi-allelic changed lines
#$multi_out - out file that lists all multi-allelic sites
#$DPR_ratio_threshold - number as proprotion e.g. 0.15 means that if 15% or more of the reads support another allele than the one that is called, then the genotype will be converted into missing data


###FORMAT field formats:
#GT:DP:ADF:ADR:AD - format of invariant site
#GT:PL:DP:ADF:ADR:AD - format of variant site


my ($input, $output, $output_changed_lines, $multi_changed_lines, $multi_out, $DPR_ratio_threshold) = @ARGV; 

open (IN, "<$input") or die "cant open $input - $!";
open (OUT, ">$output") or die "cant write to $output - $!";
open (OUT_MULTI, ">$multi_out") or die "cant write to $output - $!";
open (OUT_CH, ">$output_changed_lines") or die "can't write to $output_changed_lines - $!";
open (OUT_MULTI_CH,">$multi_changed_lines") or die "can't write to  $multi_changed_lines - $!";

my $count_new_missing = 0;
my $multiallelic_count =0;
my $multiallelic_new_missing =0;

while ( <IN> ){
    my $row = $_;
    my @output_array = ();
    if ($row =~ m/^##/) { #if row matches # in the beginning, print straight to output
      print OUT $row;
	}
    elsif ($row =~m/^#/) {
        chomp ($row);
        print OUT "$row\n";
    }
    elsif ($row !~ m/^#/) { #to ignore the lines in the beginning and only work with genotype lines
        my @cells = split;
        chomp ($row);
            if ($cells[4] =~/^\./) { #print out invariant sites
                chomp ($row);
                print OUT $row;
            }
            elsif ($cells[4]=~m/,/) { #(^A)|(^C)|(^G)|(^T)
                $multiallelic_count++;
                print OUT_MULTI "$row\n";
                #next;
                foreach my $cell (@cells){
                    if ($cell !~m/^\d:/) {
                        @output_array = (@output_array, $cell);
                    }
                    elsif ($cell =~m/^0:/){
                        my @FORMAT_array = (split(/:/,$cell));
                        if ($FORMAT_array[5]=~m/(\d+),(\d+),(\d+)/) {
                            my ($AD1,$AD2,$AD3) = ($1,$2,$3);
                            my $sum = $AD2+$AD3;
                            if ($sum < 1) { #if other alleles have no support, print out 0
                                @output_array = (@output_array, $cell);
                            }
                            elsif ($AD1 == 0) {
                                splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $multiallelic_new_missing++;
                                    print OUT_MULTI_CH "$cells[1]\t$new_cell\n";
                            }
                            else {
                                my $AD_ratio = $sum/$AD1;
                                if ($AD_ratio <=$DPR_ratio_threshold) {
                                   @output_array = (@output_array, $cell);
                                }
                                else {
                                    splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $multiallelic_new_missing++;
                                    print OUT_MULTI_CH "$cells[1]\t$new_cell\n";
                                }
                            }
                        }
                    }
                    elsif ($cell =~m/^1:/){
                        my @FORMAT_array = (split(/:/,$cell));
                        if ($FORMAT_array[5]=~m/(\d+),(\d+),(\d+)/) {
                            my ($AD1,$AD2,$AD3) = ($1,$2,$3);
                            my $sum = $AD1+$AD3;
                            if ($sum < 1) { #if other alleles have no support, print out 0
                                @output_array = (@output_array, $cell);
                            }
                            else {
                                my $AD_ratio = $sum/$AD2;
                                if ($AD_ratio <=$DPR_ratio_threshold) {
                                   @output_array = (@output_array, $cell);
                                }
                                else {
                                    splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $multiallelic_new_missing++;
                                    print OUT_MULTI_CH "$cells[1]\t$new_cell\n";
                                }
                            }
                        }
                    }
                    elsif ($cell =~m/^2:/){
                        my @FORMAT_array = (split(/:/,$cell));
                        if ($FORMAT_array[5]=~m/(\d+),(\d+),(\d+)/) {
                            my ($AD1,$AD2,$AD3) = ($1,$2,$3);
                            my $sum = $AD1+$AD2;
                            if ($sum < 1) { #if other alleles have no support, print out 0
                                @output_array = (@output_array, $cell);
                            }
                            else {
                                my $AD_ratio = $sum/$AD3;
                                if ($AD_ratio <=$DPR_ratio_threshold) {
                                   @output_array = (@output_array, $cell);
                                }
                                else {
                                    splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $multiallelic_new_missing++;
                                    print OUT_MULTI_CH "$cells[1]\t$new_cell\n";
                                }
                            }
                        }
                    }
                    else {
                        next,
                    }
                }
                } 
            else {    
                foreach my $cell (@cells){
                    if ($cell !~m/^\d:/) {
                        @output_array = (@output_array, $cell);
                    }
                    elsif ($cell =~m/^0:/) {
                        my @FORMAT_array = (split(/:/,$cell));
                        if ($FORMAT_array[5]=~m/(\d+),(\d+)/) {
                            my ($AD1,$AD2) = ($1,$2);
                            if ($AD2 == 0) {
                                @output_array = (@output_array, $cell);
                                }
                            elsif ($AD1 == 0) {
                                splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $count_new_missing++;
                                    print OUT_CH "$cells[1]\t$new_cell\n";
                            }
                            else {
                                my $AD_ratio = $AD2/$AD1;
                                if ($AD_ratio <=$DPR_ratio_threshold) {
                                    @output_array = (@output_array, $cell);
                                }
                                elsif($AD_ratio >$DPR_ratio_threshold){
                                    splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $count_new_missing++;
                                    print OUT_CH "$cells[1]\t$new_cell\n";
                                }
                                }
                            }
                            
                        }
                     elsif ($cell =~m/^1:/){
                         my @FORMAT_array = (split(/:/,$cell));
                        if ($FORMAT_array[5]=~m/(\d+),(\d+)/) {
                            my ($AD1,$AD2) = ($1,$2);
                            if ($AD1 == 0) {
                                @output_array = (@output_array, $cell);
                                }
                            elsif ($AD2 == 0) {
                                splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $count_new_missing++;
                                    print OUT_CH "$cells[1]\t$new_cell\n";
                            }
                            else {
                                my $AD_ratio = $AD1/$AD2;
                                if ($AD_ratio <=$DPR_ratio_threshold) {
                                    @output_array = (@output_array, $cell);
                                }
                                elsif($AD_ratio >$DPR_ratio_threshold){
                                    splice @FORMAT_array, 0, 1, '.';
                                    my $new_cell = join ':', @FORMAT_array;
                                    @output_array = (@output_array, $new_cell);
                                    $count_new_missing++;
                                    print OUT_CH "$cells[1]\t$new_cell\n";
                                }
                                }
                            }
                            
                        }
                     else {
                        print OUT_MULTI "wierd genotype not 0 or 1\n $row\n";
                     }
                    }
          
                }
            print OUT join ("\t",@output_array,"\n");	
       @output_array = ();
            }
    }
print "Number of new missing GTs is $count_new_missing\nNumber of new missing GTs at multiallelic sites is $multiallelic_new_missing\nNumber of multiallelic is $multiallelic_count\n";
print OUT_CH "Number of new missing GTs is $count_new_missing\nNumber of new missing GTs at multiallelic sites is $multiallelic_new_missing\nNumber of multiallelic is $multiallelic_count\n";
close IN;
close OUT;
close OUT_CH;
close OUT_MULTI_CH;
