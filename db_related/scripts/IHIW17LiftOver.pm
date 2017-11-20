package IHIW17LiftOver;
use strict;

use Exporter;

our @ISA= qw( Exporter );

our @EXPORT = qw(LiftOverLocal glTest genotypelistLiftOver);
our @EXPORT_OK = qw(LiftOverLocal glTest genotypelistLiftOver);

my $hlarefIndex=70;
my $hfilename="/home/oracle/scripts/IHIW17_AllelelistGgroups_history.txt";
open (my $hfile, "< $hfilename") or die "Could not open file '$hfilename' $!";

my $header=<$hfile>;
chop $header;
my @headers= split /\t/,$header;
my %hlahash=();
my @refAlleles=();
while (<$hfile>){
	chop;
	my @tmp=split /\t/;
	push @refAlleles,"HLA-$tmp[$hlarefIndex]";
	for(my $i=1;$i<=$#tmp;$i++){
		$hlahash{$headers[$i]}{"HLA-$tmp[$i]"}=$#refAlleles;
	}



}
sub LiftOverLocal($$$){
	my($gls,$ver,$newver)=@_;
	my @newgls=();
	for(my $i=0;$i<=$#{$gls};$i++){
		push @newgls,genotypelistLiftOver($gls->[$i],$ver);
	}
	@newgls;
}
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
sub alleleLiftOver($$){
	my($allele,$version)=@_;
	die "No such HLA version $version" unless(defined $hlahash{$version});
	die "No such allele $allele in  HLA version $version" unless(defined $hlahash{$version}{$allele});
	return $refAlleles[$hlahash{$version}{$allele}];
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
	$version=versionFormat($version);
	return $genotypelist if($version eq '3250');
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
