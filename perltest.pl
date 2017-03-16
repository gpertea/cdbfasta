#!/usr/bin/perl
use strict;
my $usage=q/.\/perltest.pl <fasta.cidx> <key>
 Will test the -P option of cdbyank to use the retrieved
 offset and show the fasta record at that position
/;
my ($cidx, $key)=@ARGV;
die($usage."\n") unless $key && -f $cidx;
my $file=$cidx;
$file=~s/\.cidx$//;
open(BIGFILE, $file)
 || die "Cannot open data file $file (for index $cidx)!\n";

my $ofs=`cdbyank -P -a '$key' $cidx`;
chomp($ofs);
die("Error: key $key not found in $cidx\n") unless length($ofs)>0;
$ofs=int($ofs);

seek(BIGFILE, 0, 2);

my $len=tell(BIGFILE);

print "--> File size is $len\n";
print "--> Record at position $ofs:\n";
seek(BIGFILE, $ofs, 0);
my $rec;
while (<BIGFILE>) {
 last if (/^>/ && $rec);
 $rec=$_;
 print $_;
}

close(BIGFILE);

