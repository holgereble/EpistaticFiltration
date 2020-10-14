use application "fan";
use application "polytope";

print STDOUT "Type 'run_order(\$dim, \$no_critical_interactions_to_be_sampled)' to run the experiment. \nType 'info' to show how the experiment works.\n";

sub run_order{
    my ($dim, $no_critical_interactions_to_be_sampled) = @_;
    my $seed = 10456; # set some seed value for reproducing the experiment
    my $path = join("", "/tmp/synthexp_dim",$dim, "_order.tex");
    my $csvpath = join("", "/tmp/synthexp_dim",$dim, "_order.csv");
    open OUT, "> $path";
    open CSV, "> $csvpath";
    print CSV join(",", "significance", "order", "eweight"), "\n";
    init_tex();
    my $no_fils = 1;
    my $no_crit_interactions=1;
    my $sd=0.1; # Set lowest standard deviation to start with. 
    my $fil;
    print STDOUT "no. critical interactions (no. filtrations): ";
    while($no_crit_interactions < $no_critical_interactions_to_be_sampled){
        ($fil,$seed)=sample_given_sd($dim, $sd, $seed);
        my $color="red";
        for (my $i=0; $i<$fil->CRITICAL_SIGNIFICANCES->size; ++$i){
            if ($fil->CRITICAL_SIGNIFICANCES->[$i] eq "significant"){
                $color="blue";
            }
            my $order = $fil->INTERACTION_ORDERS->[$i];
            my $eweight = $fil->CRITICAL_EPISTATIC_WEIGHTS->[$i]; 
            print OUT "\\draw [$color, thick] ($order , $eweight) node {\\tiny \$\\bullet\$};\n";
            print CSV join (",", $color, $order,$eweight),"\n";
            ++$no_crit_interactions;
        }
        $sd = $sd + 0.1; # Augment standard deviation for each iteration step (filtration sampling).
        print STDOUT $no_crit_interactions, " ($no_fils),", " ";
        ++$no_fils;
    }
    print CSV "\n\n\ number of filtrations generated, ", $no_fils;
    end_tex();
    print STDOUT "number of filtrations generated: ", $no_fils-1, "\n";
    print STDOUT "number of critical interactions evaluated: ", $no_crit_interactions, "\n";
    print STDOUT join("","Results stored in /tmp/synthexp_dim",$dim,"_order.tex (resp. csv).\n");
    close OUT;
    close CSV;
}

sub sample_given_sd{
    my ($dim, $sd, $seed) = @_;
    my @measurements = ();
    my $points = cube($dim)->VERTICES;
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
    return ($filtration,$seed);
}

sub init_tex{
    print  OUT "\\documentclass[crop,tikz]{standalone} 
        \\begin{document} 

    \\begin{tikzpicture}[yscale=1, scale=1, xscale=1]
        \\tikzset{
            tick/.style = {black, thick},
                tick_thin/.style = {black, thin}}
    \\draw[line width=0.2mm, gray ,step=1, dotted] (0,0) grid (6,10);
    \\foreach \\x in {0,...,6} {
        \\draw[tick,black,thin] (\\x,-0.05) -- (\\x,0); 
        \\node at (\\x,-0.1) [below,black] { \\small \$\\x\$ };}

    \\foreach \\y in {0,...,10} {
        \\draw[tick,black,thin] (0,\\y) -- (0.1,\\y);
        \\node at (-0.1,\\y) [left,black] {\\small \$\\y\$ };}
    \\node at (2,-1) [left,black] {order};
    \\node at (-1,2) [rotate=90] {epistatic weight};\n\n";
}

sub end_tex{
    print OUT "  \\end{tikzpicture}
    \\end{document}";
}

sub info_order{

print STDOUT "Arguments: \$dim, \$no_critical_interactions_to_be_sampled
Example: 'run_order(3,100)'
Description: Synthetic height functions over the \$dim-dimensional cube are generated until the number \$no_critical_interactions_to_be_sampled of critical interactions is reached. The standard deviation of the sampled height function is first set to 0.1 and augments by 0.1 each filtration sampling step. The heights of the wild type 0...0 and the standard vector e_{\$dim} are sampled with mean 53, all the other vertices with mean 50. The diagram which is drawn (TeX file) relates epistatic weights (statistical significance) to interaction order.\n"




}


