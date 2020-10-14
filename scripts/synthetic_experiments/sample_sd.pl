use application "fan";
use application "polytope";

print STDOUT "Type 'run_sd(\$dim, \$upper_sd_bound)' to run the experiment. \nType 'info' to show how the experiment works.\n";


sub run_sd{
    my ($dim, $upper_sd_bound) = @_;
    my $path = join("", "/tmp/synthexp_dim",$dim, ".tex");
    open OUT, "> $path";
    my $seed=10456;
    init_tex($upper_sd_bound);
    my $sd=0.1;
    my $fil;
    while($sd<$upper_sd_bound){
        ($fil,$seed)=sample_given_sd($dim, $sd, $seed);
        my $cnt_signif=0;
        my $color="red";
        my $crit_vertset;
        for (my $i=0; $i<$fil->CRITICAL_SIGNIFICANCES->size; ++$i){
            if ($fil->CRITICAL_SIGNIFICANCES->[$i] eq "significant"){
                ++$cnt_signif;
                my $edge_no = $fil->CRITICAL_EDGES->[$i];
                $crit_vertset = $fil->BIPYRAMIDS->[$edge_no]->first + $fil->BIPYRAMIDS->[$edge_no]->second;
            }
        }
        if($fil->CRITICAL_SIGNIFICANCES->[$fil->CRITICAL_SIGNIFICANCES->size-1] eq "significant" && ($crit_vertset->contains(0) && $crit_vertset->contains($dim+1))){$color="blue"}
        print $sd,":\n", $cnt_signif, " ", $crit_vertset," $color", "\n";
        print OUT "\\draw [$color, thick] ($sd , $cnt_signif) node {\\tiny \$\\bullet\$};\n";
        $sd = $sd + 0.1;
    }
    end_tex();
    close OUT;
    print STDOUT join("","Results stored in ", $path, ".\n")
}

sub sample_given_sd{
    my ($dim, $sd, $seed) = @_;
    my @measurements = ();
    my $points_old = load_data(join ("", "cube", $dim, ".dat"));
    my $points = gen_cube($dim);
    print "pointcheck: ", $points_old == $points, "\n";
    for (my $i=0; $i<$points->rows; ++$i){
        my $bottom = join ("", "set.seed(",$seed++,");sample = rnorm(n=100, m=50, sd=", $sd, ");sample;");
        my $top = join ("", "set.seed(",$seed++,");sample = rnorm(n=100, m=53, sd=", $sd, ");sample;");
        if ($i==0 || $i==$dim+1){ `echo "$top" > Rcommand.dat`;}
        else {`echo "$bottom" > Rcommand.dat`}

        my $Rfile = "Rcommand.dat";
        my $result = `R --no-save --quiet  < $Rfile 2>/dev/null`;
        $result =~ s/(.*\;)//g;
        $result =~ s/\n\s*//g;
        $result =~ s/(\D*\[\d*])//g;
        $result =~ s/>//g;
        my $meas = new Array<Rational>($result);
        push @measurements, $meas;
    }
    my $measurements = new Array<Array<Rational>>([@measurements]);
    my $filtration = new Filtration(POINTS=>$points, MEASUREMENTS=>$measurements);
    return ($filtration, $seed);
}


sub init_tex{
    my ($upper_sd_bound)=@_;
    print  OUT "\\documentclass[crop,tikz]{standalone} 
        \\begin{document} 


    \\begin{tikzpicture}[yscale=1, scale=1, xscale=1]
        \\tikzset{
            tick/.style = {black, thick},
                tick_thin/.style = {black, thin}}
    \\draw[line width=0.2mm, gray ,step=1, dotted] (0,0) grid ($upper_sd_bound,10);
    \\foreach \\x in {0,1,...,$upper_sd_bound} {
        \\draw[tick,black,thin] (\\x,-0.05) -- (\\x,0); 
        \\node at (\\x,-0.1) [below,black] { \\small \$\\x\$ };}

    \\foreach \\y in {0,...,10} {
        \\draw[tick,black,thin] (0,\\y) -- (0.1,\\y);
        \\node at (-0.1,\\y) [left,black] {\\small \$\\y\$ };}
    \\node at (4,-1) [left,black] {standard deviation};
    \\node at (-1,2) [rotate=90] {\\# of significant interactions};";
}

sub end_tex{

    print OUT "  \\end{tikzpicture}
    \\end{document}";
}

sub info_sd{
print STDOUT "Arguments: \$dim, \$upper_sd_bound
Example: 'run_sd(4,20)'
Description: Synthetic height functions over the \$dim-dimensional cube are generated. The standard deviation of the sampled height function is first set to 0.1 and augments by 0.1 each filtration sampling step until it reaches the value \$upper_sd_bound. The heights of the wild type 0...0 and the standard vector e_{\$dim} are sampled with mean 53, all the other vertices with mean 50. For every sampled filtration a dot is drawn in the TeX-diagram: it is blue if the special genotypes (wild type and e_{\$dim}) are seen by the last significant interaction and it is red else. The diagram relates the number of significant interactions to the standard deviations in this context.\n";
}

sub gen_cube{
    my $dim = $_[0];
    my @matrix = ();
    my $row_num=0;
    for (my $i=0; $i<=$dim; ++$i){
        foreach my $subset (@{all_subsets_of_k([(1..$dim)],$i)}){
            my @ones = (0);
            push @ones, @{$subset}; 
            my @row = (0) x ($dim+1);
            for (@ones){$row[$_]=1}
            push @{$matrix[$row_num++]}, @row;
        }
    }
    return new Matrix([@matrix]);
}

