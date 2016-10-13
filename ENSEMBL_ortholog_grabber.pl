#Written by Tristan Eversole while at the Ahituv lab at UCSF.
#Special thanks to my mentor Aaron Hardin, who provided a great deal of support and advice,
#without which this script could not have been written.

#This script accepts a list of species residing in the ENSEMBL database and finds all the one-to-one
#orthologs shared between them, as annotated by ENSEMBL. It writes these orthologs to a single file in FASTA
#format. It requires ENSEMBL's BioPerl libraries.

#I used these resources heavily in the construction of this script. Should you wish to modify it, you may
#find them useful.

#http://uswest.ensembl.org/info/docs/api/core/core_tutorial.html
#http://uswest.ensembl.org/info/docs/api/compara/compara_tutorial.html
#http://www.ensembl.org/info/docs/Doxygen/compara-api/classes.html
#http://qntm.org/files/perl/perl.html

use strict;
use warnings;
use diagnostics;
use Tie::RefHash; #This allows us to assign objects as hash keys
use List::MoreUtils; #This provides us with an easy means of removing duplicate entries from an array

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Feature;
use Bio::SeqIO;

#Set up the ENSEMBL registry so that we can do database queries
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host		=> 'ensembldb.ensembl.org',
	-user		=> 'anonymous',
	-port		=> '3306',
	-verbose	=> 1
);

Bio::EnsEMBL::Registry->set_reconnect_when_lost(1); #ENSEMBL has a bad habit of dropping connections

#This is the list of species to get orthologs from
my @common_species_names = (
	"Anolis",
	"Xenopus",
	"Chicken",
	"Zebra finch",
	"Chinese softshell turtle",
	"Human",
	"Dog",
	"Opossum",
	"Platypus",
);

#Translate between common names and whatever ENSEMBL understands by querying GenomeDB objects
my $genome_db_adaptor = $registry->get_adaptor( 'Multi', 'Compara', 'GenomeDB');
my @species_names;
my %species_name_translator;
foreach my $entry (@common_species_names) {
	my $species_genome_db = $genome_db_adaptor->fetch_by_registry_name($entry);
	my $new_species_name = $species_genome_db->name();
	#Create a hash that can translate between common names and whatever ENSEMBL understands
	$species_name_translator{$new_species_name} = $entry;
	#Create an array of names that we can use in the program
	push @species_names, $new_species_name;
};

#Use the first species in the list as our reference for performing the other queries
#If these are one-to-one orthologs, the choice of species shouldn't matter
#The second species will be used for our pairwise comparison
my $first_query_genome = $species_names[0];
my $second_query_genome = $species_names[1];
tie my %observed_orthologs, 'Tie::RefHash';
my $first_query_ortholog_count = 0;
my $second_query_ortholog_count = 0;
my $observed_ortholog_count = 0;
my %species_ortholog_counts = ();

my $method_link_species_set_adaptor = $registry->get_adaptor( 'Multi', 'Compara', 'MethodLinkSpeciesSet');
my $target_species_mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_registry_aliases("ENSEMBL_ORTHOLOGUES", [$first_query_genome,
	$second_query_genome]);

my $homology_adaptor = $registry->get_adaptor( 'Multi', 'Compara', 'Homology');
my $pairwise_orthologs = $homology_adaptor->fetch_all_by_MethodLinkSpeciesSet($target_species_mlss, "ortholog_one2one");

#Perform the first pairwise comparison

print "Fetching orthologs from $species_name_translator{$first_query_genome} and $species_name_translator{$second_query_genome}" . "..." . "\n";
# PAIRWISE: foreach my $entry (@{$pairwise_orthologs}) {
foreach my $entry (@{$pairwise_orthologs}) {
				my $gene_group = $entry->get_all_GeneMembers();
				foreach my $gene_member (@{$gene_group}) {	#There really shouldn't be more than two, ever, but you never know
					my $gene = $gene_member->get_Gene();
					#Sadly, a GeneMember object has no idea what species it's from or what type it is, requiring two extra SQL queries
					my $biotype = $gene->biotype();
					my $ortholog_species_name = $gene->species();
					if ($biotype eq 'protein_coding') {
						if ($ortholog_species_name eq $first_query_genome) {
							$observed_ortholog_count = $observed_ortholog_count + (scalar @{$gene_group});
							$first_query_ortholog_count++;
							print "Current ortholog count: $observed_ortholog_count" . "\r";
							$observed_orthologs{ $gene_member } = [ @{$gene_group } ] #This should be a reference to an array of GeneMember objects
						} elsif (($biotype eq 'protein_coding') && ($ortholog_species_name eq $second_query_genome)) {
							$second_query_ortholog_count++;							
						} else {
							print "\n" . "Error: Unknown species $ortholog_species_name found in pairwise comparison" . "\n";
						};
					};				
					# if ($observed_ortholog_count == 100) {
					# 	last PAIRWISE;
					# };
				};
		};
print "\n";

$species_ortholog_counts{$first_query_genome} = $first_query_ortholog_count;
print $species_name_translator{$first_query_genome} . " orthologs found: " . $first_query_ortholog_count . "\n";
$species_ortholog_counts{$second_query_genome} = $second_query_ortholog_count;
print $species_name_translator{$second_query_genome} . " orthologs found: " . $second_query_ortholog_count . "\n";
print "Current cumulative ortholog count: " . $observed_ortholog_count . "\n";

#Perform the other pairwise comparisons, using the first entry in @species_names as a reference

my @query_species_names = @species_names;	# @species_names will be useful later, so we don't want to destroy it now
my $reference_species = shift @query_species_names;
shift @query_species_names; #Remove the second entry, because we've already compared the first and second species

#Iterate over the array of species names

while (@query_species_names) {
	my $current_species_query = shift @query_species_names;
	print "Currently querying $species_name_translator{$reference_species} and $species_name_translator{$current_species_query}" . "..." . "\n";
	my $current_species_observed_ortholog_count = 0;
	my @reference_species_gene_list = sort keys %observed_orthologs;
	
	#Iterate over an array reference species GeneMember items
	while (@reference_species_gene_list) {
		my $reference_species_gene_member = shift @reference_species_gene_list;
		my $sequential_orthologs = $homology_adaptor->fetch_all_by_Member_paired_species($reference_species_gene_member, $current_species_query,
			['ENSEMBL_ORTHOLOGUES']);
		#If $sequential_orthologs is undefined, then this GeneMember has no homology between the two query species
		if ((defined $sequential_orthologs) && (@{$sequential_orthologs})) {
			foreach my $entry (@{$sequential_orthologs}) {
				#Retain only one-to-one orthologs-- I wish there was a way to fetch only those in the first place
				if ($entry->description() eq 'ortholog_one2one') {
					my $gene_group = $entry->get_all_GeneMembers();
					foreach my $gene_member (@{$gene_group}) {
						my $gene = $gene_member->get_Gene();
						my $biotype = $gene->biotype();
						my $ortholog_species_name = $gene->species();
						if (($biotype eq 'protein_coding') && ($ortholog_species_name eq $current_species_query)) {
							$observed_ortholog_count++;
							$current_species_observed_ortholog_count++;
							# print $ortholog_species_name . "\n";
							print $species_name_translator{$current_species_query} . " ortholog count: $current_species_observed_ortholog_count" . "\r";
							unshift( @{ $observed_orthologs{ $reference_species_gene_member } }, $gene_member);
						} elsif (($ortholog_species_name ne $reference_species) && ($ortholog_species_name ne $current_species_query)) {
							print "\n" . "Error: Unknown species $ortholog_species_name found in pairwise comparison" . "\n";
						} elsif ($biotype ne 'protein_coding') {
							$observed_ortholog_count = $observed_ortholog_count - (scalar (@{ $observed_orthologs{$reference_species_gene_member} }) );						
							delete $observed_orthologs{$reference_species_gene_member};
						};
					};
				#If it's not a one-to-one ortholog, remove it from the hash 
				#and decrement $observed_ortholog_count by the size of the array associated with the key			
				} else {
					# print "Beginning non-one-to-one deletion..." . "\n";
					# print "@{ $observed_orthologs{$reference_species_gene_member} }" . "\n";
					# print scalar @{ $observed_orthologs{$reference_species_gene_member} } . "\n";
					$observed_ortholog_count = $observed_ortholog_count - (scalar (@{ $observed_orthologs{$reference_species_gene_member} }) );						
					delete $observed_orthologs{$reference_species_gene_member};
					last;	#If we delete $observed_orthologs{$reference_species_gene_member}, we can't iterate over $sequential_orthologs anymore
				};
			};
		#If there's no homology, remove it from the hash
		#and decrement $observed_ortholog_count by the size of the array associated with the key
		} else {
			$observed_ortholog_count = $observed_ortholog_count - (scalar (@{ $observed_orthologs{$reference_species_gene_member} }) );						
			delete $observed_orthologs{$reference_species_gene_member};
		};
	};
	$species_ortholog_counts{$current_species_query} = $current_species_observed_ortholog_count;
	print $species_name_translator{$current_species_query} . " orthologs found: " . $current_species_observed_ortholog_count . "\n";
	print "Current cumulative ortholog count: " . $observed_ortholog_count . "\n";
	# my $size = scalar (keys %observed_orthologs);
	# print "Observed ortholog key number: " . $size . "\n";
};

foreach my $name (sort (keys %species_ortholog_counts)) {
	print $species_name_translator{$name} . " ortholog count: " . $species_ortholog_counts{$name} . "\n";
};

print "\n" . "Total ortholog count: " . $observed_ortholog_count . "\n";

#Finally, fetch the relevant information and write it to a file in FASTA format

print "Writing to file..." . "\n";
my $filename = 'orthologs.fa';
open (my $fh, '>', $filename) or die "Could not open file '$filename' $!"; 

my $gene_member_adaptor = $registry->get_adaptor( 'Multi', 'Compara', 'GeneMember');
my $ortholog_written_count = 0;

#Using a while-shift loop will take up less memory than a foreach loop, I think
my @final_orthologs = sort (keys %observed_orthologs);
while (@final_orthologs) {
	my $entry = shift @final_orthologs;
	my @ortholog_chain = @{ $observed_orthologs{$entry} };
	my $size = scalar @ortholog_chain;
	# print $size . "\n";
	while (@ortholog_chain) {
		my $gene_member = shift @ortholog_chain;
		my $gene = $gene_member->get_Gene();
		my $species = $gene->species();
		my $stable_id = $gene->stable_id();
		my $description = $gene_member->description();
		my $biotype = $gene->biotype();
		my $strand = $gene->strand();
		my $canonical_sequence = $gene_member->get_canonical_SeqMember();
		my $canonical_CDS = $canonical_sequence->other_sequence('cds');
		# print $species . " " . $stable_id . "\n";
		my $FASTA_string = ">" . $species . " " . $stable_id . " " . $description . " " . $biotype . " " . $strand . "\n" . $canonical_CDS . "\n";
		print $fh $FASTA_string;
		$ortholog_written_count++;
		print "Total orthologs written: " . $ortholog_written_count . "\r";
	};
};

close $fh;
print "\n" . "Done." . "\n";
