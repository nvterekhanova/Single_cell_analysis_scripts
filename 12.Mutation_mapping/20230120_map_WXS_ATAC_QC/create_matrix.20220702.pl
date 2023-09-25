#Based on create_matrix.20220626.pl
use strict;

my $sample=$ARGV[0];
my $disease=$ARGV[1];
my @all_barcodes;
my $c_str=0;
open CELL_LIST, "out/$disease/barcodes/$sample.ATAC_barcodes.tsv";
while (<CELL_LIST>){
    chomp;
    my $str=$_;
    push @all_barcodes, $_;
}

my %m_id;
my $m_str=0;
open MAF, "inputs/$disease/$disease\_PanCanAllNonSilent.maf";
while (<MAF>){
    chomp;
    my $str=$_;
    $m_str++;
    if ($m_str==1){
	next;
    }
    my @arr=split /\t/, $str;
    my $coord=$arr[4].'_'.$arr[5];
    $m_id{$coord}=$arr[0].'_'.$coord;
}


my (%barc_ref, %barc_var, @all_ids, %id_ex, $id);
open BARCODES, "out/$disease/muts.mapped/$sample/$sample.step2.out";
 LINE: while (<BARCODES>){
     chomp;
     my $str=$_;
     my @arr=split /\t/, $str;
     my ($chr, $posit);
     if ($str=~/^chr[0-9A-Z]+\s+.*Ref-Support/){
	 $chr=$arr[0];
	 $posit=$arr[1];
	 $id=$chr.'_'.$posit;
	 if ($id_ex{$id} !=1){
	     push @all_ids, $id;
	     $id_ex{$id}=1;
	 }
	 while (<BARCODES>){
	     my $ref_str=$_;
	     chomp($ref_str);
	     if ($ref_str !~/Var-Support/){
		 my @arr_ref=split /\t/, $ref_str;
 		 push @{$barc_ref{$id}}, $arr_ref[1];
	     }
	     else {
		 next LINE;
	     }
	 }
     }
     else {
	 if ($str !~/^chr/){
	     my @arr_var=split /\t/, $str;
	     push @{$barc_var{$id}}, $arr_var[1];
	     print "$id\t$arr_var[1]\n";
	 }
     }
}

open OUT, ">out/$disease/muts_assays/res_MutAssay_$sample.txt";
for (@all_barcodes){
    print OUT "\t$_";
}
print OUT "\n";
for (@all_ids){
    my $id=$_;
    print OUT "$m_id{$id}";
    for (@all_barcodes){
	my $barc=$_;
	my $found=0;
	for (@{$barc_ref{$id}}){
	    if ($_ eq $barc){
		if ($found !=1){
		    print OUT "\t1";
		    $found=1;
		}
	    }
	}
	for (@{$barc_var{$id}}){
	    if ($_ eq $barc){
		if ($found !=1){
		    print OUT "\t2";
		    $found=1;
		}
	    }
	}
	if ($found==0){
	    print OUT "\t0";
	}
    }
    print OUT "\n";
}
