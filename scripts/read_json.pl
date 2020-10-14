# polymake script to read data files from
# https://github.com/harmslab/genotype-phenotype-maps.git

use application "common";
use JSON;

# @example:
# > $j = read_json("bridgham.json");
# > $plhs=data_from_genotypes($j);
# > $sd = new fan::SubdivisionOfPoints(POINTS=>$plhs->[0], POINT_LABELS=>$plhs->[1], WEIGHTS=>$plhs->[2]);
# > print $sd->N_MAXIMAL_CELLS;
# | 19

# Input: JSON file name
# Returns: hash (reference); this is a nested hash/array type data structure.
sub read_json($) {
  my ($filename)=@_;
  open my $in, "<:utf8", $filename;
  local $/;
  my $json = JSON->new->relaxed->utf8->decode(<$in>) or die "can't read $filename: $!\n";
  close $in;
  return $json;
}

# Input: JSON hash/array reference.
# Returns: hash (reference), mapping genotypes to phenotypes.
# The genotypes are strings, where each character corresponds to one allele.
# The phenotypes are floats.
sub map_geno_pheno($) {
  my ($j)=@_;
  my $map={};
  my $genotypes=$j->{data}->{genotypes};
  my $phenotypes=$j->{data}->{phenotypes};
  my $n=scalar(@{$genotypes});
  die "genotype/phenotype mismatch" unless scalar(@{$phenotypes})==$n;
  for (my $i=0; $i<$n; ++$i) {
    $map->{$genotypes->[$i]} = $phenotypes->[$i];
  }
  return $map;
}


sub map_geno_stdev($) {
  my ($j)=@_;
  my $map={};
  my $genotypes=$j->{data}->{genotypes};
  my $stdevs=$j->{data}->{stdeviations};
  my $n=scalar(@{$genotypes});
  die "genotype/stdev mismatch" unless scalar(@{$stdevs})==$n;
  for (my $i=0; $i<$n; ++$i) {
    $map->{$genotypes->[$i]} = $stdevs->[$i];
  }
  return $map;
}

# Input: JSON hash/array reference.
# Returns: hash (reference), where keys are numbers 0,1,...,(number of alleles - 1);
# Each allel allows for an arbitrary number of mutations; the values are suitable vectors (which form the vertices of a simplex).
# No homogeneization here.
# Could, of course, also be an array; corresponding json bit is a hash though.
sub mutation_vectors($) {
  my ($j) = @_;
  my $vectors = {};
  foreach (keys %{$j->{mutations}}) {
    my $m = $j->{mutations}->{$_};
    my $n = scalar(@$m);
    my $matrix = new Array<Vector<Int>>(rows(unit_matrix<Int>($n)->minor(All,~[0])));
    my $v = {};
    for (my $i=0; $i<$n; ++$i) { 
      $v->{$m->[$i]} = $matrix->[$i];
    }
    $vectors->{$_} = $v;
  }
  return $vectors;
}

# Input: JSON hash/array reference.
# Returns: array (reference), with four entries: points, labels, height function, standard errors (suitable polymake small types).
# The points are concatenated from the mutation_vectors by reading the genotype character by character.
# Thus the points form the vertices in a product of simplices, one per allele; dimensions depend on the number of mutations.
# Homogenized coordinates.
sub data_from_genotypes($) {
  my ($j) = @_;
  my $genotypes = $j->{data}->{genotypes};
  my $n = scalar(@$genotypes);
  my $n_replicates = $j->{data}->{n_replicates};
  my $pheno_map = map_geno_pheno($j);
  my $stdev_map = map_geno_stdev($j);
  my $v = mutation_vectors($j);
  my $points = new Array<Vector>($n);
  my $labels = new Array<String>($n);
  my $heights = new Vector<Float>($n);
  my $stdevs = new Vector<Float>($n);
  for (my $i=0; $i<$n; ++$i) {
    my $g = $genotypes->[$i];
    $labels->[$i] = $g;
    my @geno = split //, $g;
    my $p = "1"; # homogenizing coordinate
    for (my $k=0; $k<scalar(@geno); ++$k) {
      $p .= " " . $v->{$k}->{$geno[$k]};
    }
    $points->[$i] = new Vector($p);
    $heights->[$i] = $pheno_map->{$g};
    $stdevs->[$i] = $stdev_map->{$g};
  } 
  return [$points, $labels, $heights, $stdevs/(sqrt($n_replicates))];
}

sub filtration_from_json{
my ($filename)=@_;
my $j=read_json($filename);
my $plhs=data_from_genotypes($j);
my $points_matrix=new Matrix<Float>($plhs->[0]);
my $setope = $points_matrix|$plhs->[3];
use application "polytope";
$filename =~ s/.json//g;
my $fil = new Filtration($filename, POINTS=>$points_matrix, POINT_LABELS=>$plhs->[1], HEIGHTS=>-$plhs->[2], SETOPE=>$setope);
use application "common";
return $fil;
}
