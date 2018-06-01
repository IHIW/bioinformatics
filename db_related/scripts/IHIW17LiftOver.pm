package IHIW17LiftOver;
use strict;

use Exporter;

our @ISA= qw( Exporter );

our @EXPORT = qw(LiftOverLocal glTest glLiftOver genotypelistLiftOver);
our @EXPORT_OK = qw(LiftOverLocal glTest glLiftOver genotypelistLiftOver);

my $hlarefIndex=70;
my $hfilename="/home/oracle/scripts/IHIW17_AllelelistGgroups_history.txt";
open (my $hfile, "< $hfilename") or die "Could not open file '$hfilename' $!";

my $header=<$hfile>;
chop $header;
my @headers= split /\t/,$header;
my %hlahash=();
my %digit6hash=();
my %digit4hash=();
my @refAlleles=();
my %refAlleleIndex;
while (<$hfile>){
	chop;
	my @tmp=split /\t/;
	push @refAlleles,"HLA-$tmp[$hlarefIndex]";
	$refAlleleIndex{"HLA-$tmp[$hlarefIndex]"}=$#refAlleles unless(defined $refAlleleIndex{"HLA-$tmp[$hlarefIndex]"});
	for(my $i=1;$i<=$#tmp;$i++){
		if($tmp[$i]=~/^(.*:.*:.*):.*[^G]/){
			my $coding=$1;
			$coding=~/^(.*:.*):.*/;
			my $digit4=$1;
			if(defined $digit6hash{$headers[$i]}{"HLA-$coding"}){
				 $digit6hash{$headers[$i]}{"HLA-$coding"}.="/HLA-$tmp[$i]" unless(defined $hlahash{$headers[$i]}{"HLA-$tmp[$i]"});
			}
			else{
				$digit6hash{$headers[$i]}{"HLA-$coding"}="HLA-$tmp[$i]";
				
			}
	
			if(defined $digit4hash{$headers[$i]}{"HLA-$digit4"}){
				 $digit4hash{$headers[$i]}{"HLA-$digit4"}.="/HLA-$tmp[$i]" unless(defined $hlahash{$headers[$i]}{"HLA-$tmp[$i]"});
			}
			else{
				$digit4hash{$headers[$i]}{"HLA-$digit4"}="HLA-$tmp[$i]";
				
			}
		}
		$hlahash{$headers[$i]}{"HLA-$tmp[$i]"}=$#refAlleles;
	}



}
sub ExpendGL($$){
	my($glstring,$version)=@_;
	my @alleles=split /[\|\^\+\/\~]/,$glstring;
	my @seperator=split /[^\|\^\+\/\~]+/,$glstring;
	my $newglstring="";
	for(my $i=0;$i<=$#alleles;$i++){
		if(defined $hlahash{$version}{$alleles[$i]}){
			$newglstring.=$alleles[$i];
		}
		else{
			if(defined $digit6hash{$version}{$alleles[$i]}){
				$newglstring.=$digit6hash{$version}{$alleles[$i]};
			}
			else {
				if(defined $digit4hash{$version}{$alleles[$i]}){
					$newglstring.=$digit4hash{$version}{$alleles[$i]};
				}
				else{
				die "No such allele nor corresponding 8-digits alleles for $alleles[$i] in  HLA version $version";
				}
			}
		}
		$newglstring.=$seperator[$i+1];
	}
	return $newglstring;
}
sub LiftOverLocal($$$){
	my($gls,$ver,$newver)=@_;
	my @newgls=();
	for(my $i=0;$i<=$#{$gls};$i++){
		push @newgls,glLiftOver($gls->[$i],$ver);
	}
	@newgls;
}
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
sub alleleLiftOver($$){
	my($number,$version)=@_;
	return $refAlleleIndex{$refAlleles[$number]};
}
sub allelelistLiftOver($$){
	my($allelelist,$version)=@_;
	my @alist=split /\//,$allelelist;
	my @liftlist=();
	push @liftlist,alleleLiftOver($_,$version) foreach(@alist);
	my @matchlist=();
	for(my $i=0;$i<=$#alist;$i++){
		push @matchlist,$alist[$i] if($alist[$i] eq $liftlist[$i]);
	}
	if($#matchlist>=0){
		return join('/',@matchlist);
	}
	else{
		
		return join('/',uniq(@liftlist));
	}

}
sub genotypeLiftOver($$){
	my($genotype,$version)=@_;
	my @liftlist=();
	push @liftlist, allelelistLiftOver($_,$version) foreach(split /\+/,$genotype);
	return join('+',@liftlist);

}
sub genotypelistLiftOver($$){
	my($genotypelist,$version)=@_;
	my @gllist=split /\|/,$genotypelist;
	my @liftlist=();
	push @liftlist, genotypeLiftOver($_,$version) foreach(@gllist);
	my @matchlist=();
	my @onematchlist=();
	for(my $i=0;$i<=$#gllist;$i++){
		push @matchlist,$gllist[$i] if($gllist[$i] eq $liftlist[$i]);
		push @onematchlist,$liftlist[$i] if(onematch($gllist[$i] ,$liftlist[$i]));
	}
	if($#matchlist>=0){
		push @matchlist, $_ foreach(@onematchlist);
		return join('|',uniq(@matchlist));
	}
	else{
		return join('|',uniq(@liftlist));
	}
}
sub glLiftOver($$){
	my($glstring,$version)=@_;
	$version=versionFormat($version);
	$glstring=ExpendGL($glstring,$version);
	if($version eq '3250'){
		my @alleles=split /[\|\^\+\/\~]/,$glstring;
		foreach my $allele (@alleles){
			die "No such allele $allele in  HLA version $version" unless(defined $hlahash{$version}{$allele});
		}
		return $glstring;
	}
	return "" if($glstring eq "");
	my @liftlist=();
	my $glnumber=fromGl2Number($glstring,$version);
#	print $glnumber."\n";;
	push @liftlist, genotypelistLiftOver($_,$version) foreach(split /\^/,$glnumber);
#	print join('^',@liftlist)."\n";
	return fromNumber2GL(join('^',@liftlist));
}

sub fromGl2Number($$){
	my($glstring,$version)=@_;
	my @alleles=split /[\|\^\+\/\~]/,$glstring;
	my @delimiters= ( $glstring =~ /[\|\^\+\/\~]/g );
	my @numbers=();
	push @numbers, allele2Number($_,$version) foreach(@alleles);
	my $numberstr=$numbers[0];
	for(my $i=0;$i<=$#delimiters;$i++){
		$numberstr.=$delimiters[$i].$numbers[$i+1];
	}
	return $numberstr;
	


}
sub fromNumber2GL($){
	my($glnumber,$version)=@_;
	my @numbers=split /[\|\^\+\/\~]/,$glnumber;
	my @delimiters= ( $glnumber=~ /[\|\^\+\/\~]/g );
	my @alleles=();
	push @alleles, $refAlleles[$_] foreach(@numbers);
	my $glstring=$alleles[0];
	for(my $i=0;$i<=$#delimiters;$i++){
		$glstring.=$delimiters[$i].$alleles[$i+1];
	}
	return $glstring;


}
sub allele2Number($$){
	my($allele,$version)=@_;
	die "No such HLA version $version" unless(defined $hlahash{$version});
	die "No such allele $allele in  HLA version $version" unless(defined $hlahash{$version}{$allele});
	return $hlahash{$version}{$allele};
}

sub versionFormat($){
	my($oldv)=@_;
	my $newv="3.25.0";
	if($oldv=~/(\d+\.\d+)(\.\d+)/){
		$newv="$1.0";
	}
	elsif($oldv=~/(\d+\.\d+)/){
		$newv=$1.".0";
	}
	elsif($oldv=~/(\d)(\d)(\d)(\d)/){
		$newv="$1.$2$3.0";
	}
	$newv=~s/\.//g;
	$newv;
}
sub glTest($$){
	my @glerrors=("duplicate locus block",'duplicate genotype lists','duplicate genotypes','duplicate allele lists');
	my ($glstring,$sampleid) = @_;
#	print "$sampleid $glstring\n";
	return if($glstring eq "");
	my $glcom="/home/oracle/scripts//checkglAsFun.py -g \"".$glstring."\"";
	my $glreturn = `$glcom`;
	$glreturn=~s/\n//;
	my @glreturns=split /,/,$glreturn;
	my $errmes="";
	for(my $i=0;$i<4;$i++){
		if($glreturns[$i] ne "OK"){
			$errmes=$errmes.$glerrors[$i].":".$glreturns[$i]."\n";
		}
	}
	die "sample:$sampleid glstring errors: $errmes" if($errmes ne "");
}
sub onematch($$){
	my($gllist1,$gllist2)=@_;
	my @gllist1s=split /\+/,$gllist1;
	my @gllist2s=split /\+/,$gllist2;
	if((($gllist1s[0] eq $gllist2s[0]) + ($gllist1s[1] eq $gllist2s[1]))==1){
		return 1;
	}
	else{
		return 0;
	}
}
1;
