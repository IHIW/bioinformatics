#!/usr/bin/perl -w
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

sub getSampleContext($){
	my($filename)=@_;
	open (my $ff, "< $filename") or return "";
	my $context="";
	my $tag=0;
	my $tmp="";
	while(<$ff>){
		$tag=1 if(/<sample/);
		if(/<\/sample/){
			$context.=$tmp.$_;
			$tmp="";
			$tag=2;
		}
		if($tag==1){
			$tmp.=$_;
		}
		
	}
	close $ff;
	$context;

}

die("wrong number of parameters")  if($#ARGV<1);
my $labcode=$ARGV[0];
my $folder="/w_u/$labcode/upload/".$ARGV[1];
$ARGV[1]=~/\/?([^\/]*)\/?$/;
my $newfile="/w_u/$labcode/upload/$1.concat.hml";
opendir(my $dh, $folder) || die "Can't open $folder: $!";
my @files = grep { /\.([Hh][Mm][Ll])|([Xx][Mm][Ll])/i && -f "$folder/$_"} readdir($dh);
closedir $dh;
exit if($#files<0);
my $file = "$folder/$files[0]";
open FILE,"< $file" or die "cannot open $file";
my $head;
my $end;
my $samples;
my $tag=0;
while(<FILE>){
	$tag=1 if(/<sample/);
	if(/<\/sample/){
		$tag=2;
		next;
	}
	if($tag==0){
		$head.=$_;
	}
	if($tag==2){
		$end.=$_;
	}
	
}
close FILE;
open FILE, "> $newfile" or die "cannot open $newfile: $!";
print FILE $head;
foreach(@files){
	$file = "$folder/$_";

	print FILE getSampleContext($file);
}
print FILE $end;
