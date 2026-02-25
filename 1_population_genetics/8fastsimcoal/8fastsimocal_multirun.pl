
foreach my $i (1 .. 100){
	my $rd = int(rand(100000));
	my $line = "fsc28 -t 61ind.tpl -e 61ind.est -n100000 --foldedSFS -m -M -c4 -B4 -L50 -q -x -r $rd";
        `mkdir -p 61ind$i`;
	`cp 61ind.tpl 61ind$i`;
	`cp 61ind.est 61ind$i`;
	`cp 61ind_jointMAFpop1_0.obs 61ind$i`;
        open OUT,">./61ind$i/run.sh" || die "$!\n";
        print OUT "$line\n";
        close OUT;
}
