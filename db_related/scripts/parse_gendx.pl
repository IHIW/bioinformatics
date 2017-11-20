#!/usr/bin/perl -w
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use lib '/home/oracle/scripts';
use IHIW17LiftOver;


my %alleleNames325;

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
my %uniqlc;
my $labcode=$ARGV[0];
my $inXml=$ARGV[1];
#$ARGV[1]=~/\/?([^\/]*)\/?$/;
my $reportId=$ARGV[1];
$reportId=~s/^.*\///;
$reportId=~s/\.tarr.*$//;
my $outXml="/w_u/$labcode/upload/ws_xml/$reportId.ihiw.xml";
-f $inXml or die "file $inXml doesn't exist: $!";
my $hax = XMLin($inXml,ForceArray => ['Mismatch','MatchID','Sample','Haplotypes','Haplotype','HaplotypeID','allele','Genotype','Match','Matchset','Locus','File']) or die "cannot parse $inXml: $!";
my $hao={};
$hao->{Software_Applied}={};
$hao->{Software_Applied}->{Software_Manufacturer}="GenDx";
$hao->{Software_Applied}->{Software_Functions}="Allele Calling:Base Calling/Consensus Generation:Read Processing";
$hao->{Software_Applied}->{Software_Name}=$hax->{AnalysisSoftware}->{Software};
$hao->{Software_Applied}->{Software_Version}=$hax->{AnalysisSoftware}->{Version};

my $soft_version=$hao->{Software_Applied}->{Software_Version};
$soft_version=~/^(\d+\.\d+)/;
$soft_version=$1;
my $arrsampleo=[];
my $arrsample=$hax->{Samples}->{Sample};
my $sample;
my $i=0;
foreach $sample (@{$arrsample}){
	my $has={};
	my $arrLocio=[];
	my $arrFeatureo=[];
	$has->{SampleID}=$sample->{Name};
	my $files = $sample->{Files}->{File};
	my $fileLoc = join ",", @{$files};
	my $arrLoci = $sample->{Loci}->{Locus};
	my $locus;
	my $j=0;
	my $k=0;
	my @arrgls2;
	my $dbn="3.25.0";
	foreach $locus (@{$arrLoci}){
		my $dp = $locus->{ReadDepth}->{Mean};
		my $lname = $locus->{Name};
		my $db = $locus->{AlleleDB}->{Version};
		$db=~/ (\d\.\d+)/;
		$dbn=$1.".0";
#		die("do not support version higher than 3.25")  if($dbn=~/3.2[6789]/);		
		
		#feature Locations
		my $region;
		${$arrFeatureo}[$j]={};
		my $featurenum=0;
		my @featureStarts=();
		foreach $region (@{$locus->{AlleleDB}->{Regions}->{Region}}){
			${$arrFeatureo}[$j]->{Locus} = $lname;
			${$arrFeatureo}[$j]->{FeatureStart} = $region->{begin};
			${$arrFeatureo}[$j]->{FeatureStop} = $region->{end};
			${$arrFeatureo}[$j]->{FeatureNumber} = $featurenum;
			$featurenum = $featurenum+1;
			push @featureStarts,$region->{begin};
			$j=$j+1;
		}
		#allele Info
		my %hacons=();
		my $haplo;
		if($soft_version<=2.2){
		foreach $haplo (@{$locus->{Haplotypes}->[0]->{Haplotypes}}){
			$hacons{$haplo->{ID}}={};
			$hacons{$haplo->{ID}}->{cons}=$haplo->{content};
			$hacons{$haplo->{ID}}->{begin}=$haplo->{begin};
			$hacons{$haplo->{ID}}->{end}=$haplo->{end};
			$haplo->{beginGap}=0 unless(defined $haplo->{beginGap});
			$hacons{$haplo->{ID}}->{beginGap}=$haplo->{beginGap};
			$haplo->{endGap}=0 unless(defined $haplo->{endGap});
			$hacons{$haplo->{ID}}->{endGap}=$haplo->{endGap};
		}
		}
		else{
		foreach $haplo (@{$locus->{Haplotypes}->[0]->{Haplotype}}){
			$hacons{$haplo->{ID}}={};
			$hacons{$haplo->{ID}}->{cons}=$haplo->{content};
			$hacons{$haplo->{ID}}->{begin}=$haplo->{begin};
			$hacons{$haplo->{ID}}->{end}=$haplo->{end};
			$haplo->{beginGap}=0 unless(defined $haplo->{beginGap});
			$hacons{$haplo->{ID}}->{beginGap}=$haplo->{beginGap};
			$haplo->{endGap}=0 unless(defined $haplo->{endGap});
			$hacons{$haplo->{ID}}->{endGap}=$haplo->{endGap};
		}
		}
		#mismatch Info
		my $mcb;
		my %misMatches=();
		foreach my $mset (@{$locus->{Matching}->{Matchset}}){
			$mcb = $mset->{MatchCombination};
			foreach my $matchid (@{$mcb->{MatchID}}){
				$misMatches{$matchid->{RefAllele}}=$matchid->{Mismatch};
			}
		}
		$mcb=$locus->{Matching}->{Matchset}->[0]->{MatchCombination};
		#summary
		my $match;
		foreach $match (@{$locus->{Matching}->{Matches}->{Match}}){
			my $allele=$match->{ID};
			my @phasings = split /,/,$match->{phasing};
			$phasings[0] = 1 unless($match->{phasing});
			my $hid;
			my @hapcom=();
			foreach $hid (@{$match->{HaplotypeCombination}->{HaplotypeID}}){
				my $hpkey=$hid;
				my $phasing=0;
				if($hpkey=~/heterozygous(\d+)/){
					$phasing=$phasings[$1-1];
				}
				my $uk=$sample->{Name}."_".$allele."_".$phasings[0]."_".$hpkey;
				next if(defined $uniqlc{$uk});
				$uniqlc{$uk}=1;
				push @hapcom,$hacons{$hpkey};
			
=old
				${$arrLocio}[$k]={};
				${$arrLocio}[$k]->{HLATyping}=$allele;
				${$arrLocio}[$k]->{Alignment_Reference_DB}=$db;
				${$arrLocio}[$k]->{BaseCalling_Reference_DB}=$db;
				${$arrLocio}[$k]->{Consensus_Sequence}=$hacons{$hpkey}->{cons};
				${$arrLocio}[$k]->{Start_Position}=$hacons{$hpkey}->{begin};
				${$arrLocio}[$k]->{Feature}="Genomic - Unknown Location";
				${$arrLocio}[$k]->{Locus_name}=$lname;
				${$arrLocio}[$k]->{MeanReadDepth}=$dp;
				${$arrLocio}[$k]->{DataFileLoc}=$fileLoc;
				${$arrLocio}[$k]->{PhasingGroup}=$phasing;
				${$arrLocio}[$k]->{FeatureNumber}=buildStr(findLoc($hacons{$hpkey}->{begin},\@featureStarts),findLoc($hacons{$hpkey}->{end},\@featureStarts));
=cut
			}
			my @sorted =  sort { $a->{begin} <=> $b->{begin} || $a->{end} <=> $b->{end} } @hapcom;
			my @newhapcom=({});
			my $j=0;
			foreach $hid (@sorted){
			#	$hid->{endGap}=0 unless(defined $hid->{endGap});
			#	$hid->{beginGap}=0 unless(defined $hid->{beginGap});
				if(!defined($newhapcom[$j]{end}) || ($hid->{begin}+$hid->{beginGap})==$newhapcom[$j]{end}+1+$newhapcom[$j]{endGap}){
					$newhapcom[$j]{begin} = $hid->{begin} unless(defined $newhapcom[$j]{begin});
					$newhapcom[$j]{end}=$hid->{end};
					$newhapcom[$j]{endGap}=$hid->{endGap};
					$newhapcom[$j]{cons}.=$hid->{cons};

					
				}
				else{
					$j=$j+1;
					$newhapcom[$j]={begin=>$hid->{begin},end=>$hid->{end},endGap=>$hid->{endGap},cons=>$hid->{cons}};
				}
		#		print $has->{SampleID}." ".$allele." ".$phasings[0]." ". $newhapcom[$j]{begin}." ".$newhapcom[$j]{end}." " .$newhapcom[$j]{endGap}."\n";

			}
			foreach $hid (@newhapcom){
				${$arrLocio}[$k]={};

				my $orgallele="HLA-$allele";
				if(defined $alleleNames325{$orgallele}){
					${$arrLocio}[$k]->{HLATyping}=$alleleNames325{$orgallele};
				}
				else{
					$alleleNames325{$orgallele}=&genotypelistLiftOver($orgallele,$dbn);
					${$arrLocio}[$k]->{HLATyping}=$alleleNames325{$orgallele};
				}
				my $polyStr=$db." ". $alleleNames325{$orgallele}." :";
				${$arrLocio}[$k]->{Alignment_Reference_DB}="IMGT/HLA 3.25.0";
				${$arrLocio}[$k]->{BaseCalling_Reference_DB}="IMGT/HLA 3.25.0";
				my $tmpcon=$hid->{cons};
				$tmpcon=~s/-//g;
				${$arrLocio}[$k]->{Consensus_Sequence}=$tmpcon;
				${$arrLocio}[$k]->{Start_Position}=$hid->{begin};
				${$arrLocio}[$k]->{Feature}="Genomic - Unknown Location";
				${$arrLocio}[$k]->{Locus_name}=$lname;
				${$arrLocio}[$k]->{MeanReadDepth}=$dp;
				${$arrLocio}[$k]->{DataFileLoc}=$fileLoc;
				${$arrLocio}[$k]->{PhasingGroup}=$phasings[0];
#				print "$allele ".$hid->{begin}." ".$phasings[0]."\n";
		
				if($#{$misMatches{$allele}}>=0){
					my $mismatch;
					foreach $mismatch (@{$misMatches{$allele}}){
						$mismatch->{IMGT}=~s/:\d+//;
						if($mismatch->{IMGT}>=$hid->{begin} and $mismatch->{IMGT}<=$hid->{end}){
								
							$polyStr=$polyStr.XMLout($mismatch,RootName=>"Mismatch")." ";
						}
					}
				}
				$polyStr=~s/\n//g;
				${$arrLocio}[$k]->{NovelPolymorphism}=$polyStr;
				$k=$k+1;
			}
		}
		#build GL String;
		my $genotype;
		my @arrAlleles1=();
		my @arrAlleles2=();
		unless (defined $locus->{TypingResult}->{GenoTypeList}->{Genotype}){
			my @arrAlleles= keys %misMatches;
			if($#arrAlleles==-1){
			}
			elsif($#arrAlleles==0){
				push @arrAlleles1,{$arrAlleles[0]=>1};
				push @arrAlleles2,{$arrAlleles[0]=>1};
			}
			else{
				push @arrAlleles1,{$arrAlleles[0]=>1};
				push @arrAlleles2,{$arrAlleles[1]=>1};
			}
		}
		foreach $genotype (@{$locus->{TypingResult}->{GenoTypeList}->{Genotype}}){
			my $allele;
			my ($m1,$m2);
			foreach $allele (@{$genotype->{allele}}){
				$m1+=$allele->{PriorityMM};
				$m2+=$allele->{NonPriorityMM};
			}
			if($m1==$mcb->{PriorityMM} and $m2==$mcb->{NonPriorityMM}){
				my $newPair=1;
				
				my $allele1=$genotype->{allele}->[0]->{content};
				my $allele2=$genotype->{allele}->[1]->{content};
				my $ind;
				for($ind=0;$ind<=$#arrAlleles1;$ind++){
					if(defined $arrAlleles1[$ind]{$allele1}){
						$arrAlleles2[$ind]{$allele2}=1;
						$newPair = 0;
						last;
					}
					elsif(defined $arrAlleles1[$ind]{$allele2}){
						$arrAlleles2[$ind]{$allele1}=1;
						$newPair = 0;
						last;
					}
					elsif(defined $arrAlleles2[$ind]{$allele1}){
						$arrAlleles1[$ind]{$allele2}=1;
						$newPair = 0;
						last;
					}
					elsif(defined $arrAlleles2[$ind]{$allele2}){
						$arrAlleles1[$ind]{$allele1}=1;
						$newPair = 0;
						last;
					}
					
				}
				if($newPair==1){
					$arrAlleles1[$ind]={$allele1=>1};
					$arrAlleles2[$ind]={$allele2=>1};

				}
			}
			
		}
		my @arrgls;
		for(my $ind=0;$ind<=$#arrAlleles1;$ind++){
			$arrgls[$ind]=(join "/",map { "HLA-" . $_ } keys %{$arrAlleles1[$ind]})."+".(join "/", map { "HLA-" . $_ } keys %{$arrAlleles2[$ind]});
		}
		push @arrgls2, ( join "|",@arrgls) if($#arrgls>=0);
	}
	$has->{Genotyping}={};
	$has->{Genotyping}->{Genotype_GL}=join "^", &LiftOverLocal(\@arrgls2,$dbn,"3.25.0");
	my $oldglstring= join "^", @arrgls2;
	$has->{Genotyping}->{Original_GL}=$dbn.":".$oldglstring unless($oldglstring eq $has->{Genotyping}->{Genotype_GL});
	glTest($has->{Genotyping}->{Genotype_GL},$has->{SampleID});
	$has->{Genotyping}->{Locus}=$arrLocio;
#	$has->{Genotyping}->{FeatureCoordinate}=$arrFeatureo;
	${$arrsampleo}[$i]=$has;
	$i = $i+1;
}
$hao->{Sample}=$arrsampleo;
my $finalha={};
$finalha->{Lab}=$hao;
$finalha->{Lab}->{LabCode}=$labcode;
$finalha->{Lab}->{Lab_defined_ID}=$reportId;
open my $fhxml, '>:encoding(utf8)', $outXml or die "cannot open ($outXml): $!";
XMLout($finalha,NoSort=>1,RootName=>"IHIW_Report",OutputFile=>$fhxml);
close $fhxml;
