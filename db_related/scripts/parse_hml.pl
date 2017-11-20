#!/usr/bin/perl
use strict;
use XML::Simple;
use Data::Dumper;
use lib '/home/oracle/scripts';
use IHIW17LiftOver;


sub lcp
{
    my $match = shift;
    substr($match, (($match ^ $_) =~ /^\0*/, $+[0])) = '' for @_;
    $match;
}
sub findLoc($$){
	my($pos,$arrref)=@_;
	for(my $i=0;$i<=$#{$arrref};$i++){
		return $i-1 if($pos<${$arrref}[$i]);
	}
	return $#{$arrref};
}
sub buildStr($$){
	my($s,$e)=@_;
	my $ss="";
	for(my $i=$s;$i<=$e;$i++){
		$ss="$ss$i+";
	}
	$ss=~ s/\+$//;
	return $ss;
}

die("wrong number of parameters")  if($#ARGV<1);
my $labcode=$ARGV[0];
my $inXml=$ARGV[1];
my $reportId=$ARGV[1];
$reportId=~s/^.*\///;
$reportId=~s/\..*$//;
my $outXml="/w_u/$labcode/upload/ws_xml/$reportId.ihiw.xml";
-f $inXml or die "file $inXml doesn't exist: $!";
#my $hax = XMLin($inXml,KeyAttr => {'reference-sequence'=>'id'},ForceArray => ['sample','typing','haploid','consensus-sequence','consensus-sequence-block','variant','allele-assignment','typing-method','reference-database','reference-sequence']) or die "cannot parse $inXml: $!";
my $hax = eval {XMLin($inXml,KeyAttr => {'reference-sequence'=>'id'},ForceArray => ['sample','typing','haploid','consensus-sequence','consensus-sequence-block','variant','allele-assignment','typing-method','reference-database','reference-sequence'])};
die "cannot parse $inXml: $@" if($@);
my $hao={};
if(defined $hax->{"project-name"} && $hax->{"project-name"} ne ""){
	$reportId=$hax->{"project-name"};
}
elsif(defined $hax->{hmlid}->{root} && $hax->{hmlid}->{root} ne ""){
	$reportId=$hax->{hmlid}->{root} unless($hax->{hmlid}->{root}=~/^\d/);
}
my $clabcode = $hax->{"reporting-center"}->{"reporting-center-id"};
#if($clabcode eq ""){
	$clabcode=$labcode;
#}
my $arrsampleo=[];
my $arrsample=$hax->{sample};
my $sample;
my $i=0;
foreach $sample  (@{$arrsample}){
	my $has={};
	my $arrLocio=[];
	$has->{SampleID}=$sample->{id};
	my $arrLoci = $sample->{typing};
	my $locus;
	my $k=0;
	my (@oldgls,@vers);
	foreach $locus (@{$arrLoci}){
		next unless defined $locus->{"allele-assignment"};
		my %halocus;
		my $db="";
		foreach my $aa (@{$locus->{"allele-assignment"}}){
			$aa->{glstring}=~s/[\n\s]//g;
			$aa->{glstring}=~s/[viex]\d+//g;
#			print $aa->{glstring}."\n";
			push @oldgls,$aa->{glstring};
			glTest($aa->{glstring},$has->{SampleID});
			push @vers,$aa->{"allele-version"};
		#	die("do not support version higher than 3.25")  if($aa->{"allele-version"}=~/3.2[6789]/);
		}
	}
#	print $#oldgls."\n";
	my @newgls=&LiftOverLocal(\@oldgls,$vers[0],"3.25.0") if($#oldgls>=0);
	foreach $locus (@{$arrLoci}){
		next unless defined $locus->{"allele-assignment"};
		my %halocus;
		my $db="";
		foreach my $aa (@{$locus->{"allele-assignment"}}){
			$aa->{glstring}=shift @newgls;
		#	print $aa->{glstring}."\n";
			$aa->{"newallele-version"}="3.25.0";
		}
	}


	my @arrgls;
	foreach $locus (@{$arrLoci}){
		next unless defined $locus->{"allele-assignment"};
		my %halocus;
		my $db="";
		my $orgdb;
		foreach my $aa (@{$locus->{"allele-assignment"}}){
		$db = $aa->{"allele-db"}." ". $aa->{"newallele-version"};
		$orgdb = $aa->{"allele-db"}." ". $aa->{"allele-version"};
		my $gls=$aa->{glstring};
		$gls=~s/[\n\s]//g;
		push @arrgls,$gls if($gls=~/^HLA/);
		foreach my $haploid (@{$aa->{haploid}}){
			my $allele=$haploid->{locus}."*".$haploid->{type};
			$halocus{$allele}=$haploid->{locus};
		}
		}
		my @pgla=split /[\+\|\^\/]/,$arrgls[$#arrgls];
		my $fileLoc = $locus->{"typing-method"}->[0]->{"sbt-ngs"}->{"raw-reads"}->{uri};
		
		#summary
		foreach my $conseq (@{$locus->{"consensus-sequence"}}){
			my $haref = {};
			foreach my $refdatabase (@{$conseq->{"reference-database"}}){
				$haref = {%$haref,%{$refdatabase->{"reference-sequence"}}} if(defined $refdatabase->{"reference-sequence"});

			}
#			my $haref = $conseq->{"reference-database"}->{"reference-sequence"};
		foreach my $block (@{$conseq->{"consensus-sequence-block"}}){
			my $sid=$block->{"reference-sequence-id"};
			my $allele=$haref->{$sid}->{name};
			$allele = $haref->{name} unless defined $allele;
			my $orgallele=$allele;
			if($#pgla==0){
				$allele=$pgla[0];
			}
			elsif($#pgla>0){
				my $maxprefix=length lcp($pgla[0],$allele);
				my $maxi=0;
				for(my $pglai=1;$pglai<=$#pgla;$pglai++){
					my $prefixi=length lcp($pgla[$pglai],$allele);
					if($prefixi>$maxprefix){
						$maxprefix=$prefixi;
						$maxi=$pglai;
					}
				}
				$allele=$pgla[$maxi];
			}
			$allele="unknown" if($allele eq "");
			
			my $variants="$orgdb $orgallele:";
			#variants
			foreach my $va (@{$block->{variant}}){
				$variants.="(reference-bases=\"".$va->{"reference-bases"};
				$variants.="\" alternate-bases=\"".$va->{"alternate-bases"};
				$variants.="\" start=\"".$va->{"start"};
				$variants.="\" end=\"".$va->{"end"}."\")";

			}
			${$arrLocio}[$k]={};
			${$arrLocio}[$k]->{HLATyping}=$allele;
			${$arrLocio}[$k]->{Alignment_Reference_DB}=$db;
			${$arrLocio}[$k]->{AlleleCalling_Reference_DB}=$db;
			my $cons=$block->{sequence};
			$cons=~s/[\n\s]//g;
			${$arrLocio}[$k]->{Consensus_Sequence}=$cons;
			${$arrLocio}[$k]->{Start_Position}=$block->{start};
			${$arrLocio}[$k]->{Feature}="Genomic - Unknown Location";
			if(defined $halocus{$allele}){
				${$arrLocio}[$k]->{Locus_name}=$halocus{$allele};
			}
			else{
				$allele=~/(.*)\*/;
				${$arrLocio}[$k]->{Locus_name}=$1;
			}
			unless ($variants eq ""){
				${$arrLocio}[$k]->{NovelPolymorphism}=$variants;
			} 
			${$arrLocio}[$k]->{DataFileLoc}=$fileLoc;
			${$arrLocio}[$k]->{PhasingGroup}=$block->{"phase-set"};
				$k=$k+1;
		}
		}
	}
	$has->{Genotyping}={};
	$has->{Genotyping}->{Genotype_GL}=join "^",@arrgls;
	my $oldglstring= join "^", @oldgls;
	$has->{Genotyping}->{Original_GL}=$vers[0].":".$oldglstring unless($oldglstring eq $has->{Genotyping}->{Genotype_GL});
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\/[^\|\+\~]*)?\/\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\/[^\|\+\~]*)?\/\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\/[^\|\+\~]*)?\/\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\/[^\|\+\~]*)?\/\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\/[^\|\+\~]*)?\/\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\|[^\+\~]*)?\|\1/$1$2/g;
#	$has->{Genotyping}->{Genotype_GL}=~s/(HLA-.*)(\|[^\+\~]*)?\|\1/$1$2/g;
#	glTest($has->{Genotyping}->{Genotype_GL},$has->{SampleID});
	$has->{Genotyping}->{Locus}=$arrLocio;
	${$arrsampleo}[$i]=$has;
	$i = $i+1;
}
$hao->{Sample}=$arrsampleo;
my $finalha={};
$finalha->{Lab}=$hao;
$finalha->{Lab}->{LabCode}=$clabcode;
$finalha->{Lab}->{Lab_defined_ID}=$reportId;
open my $fhxml, '>:encoding(utf8)', $outXml or die "cannot open ($outXml): $!";
XMLout($finalha,NoSort=>1,RootName=>"IHIW_Report",OutputFile=>$fhxml);
close $fhxml;
