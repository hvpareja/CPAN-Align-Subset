package Bio::Align::Subset;


use 5.006;
use strict;
no strict "refs";
use warnings;
use Carp;
use Bio::SeqIO;


=head1 NAME

Bio::Align::Subset - A BioPerl module to generate new alignments as subset from larger alignments

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Given an array of codon positions and an alignment, the function alignSubset
returns a new alignment with the codons at those positions from the original
alignment.

Usage example:
    
    use strict;
    use warnings;
    use Data::Dumper;
    
    use Bio::Align::Subset;
    
    # The alignment in a file
    my $filename = "alignmentfile.fas";
    # The format
    my $format = "fasta";
    
    # The subset of codons
    my $subset = [1,12,25,34,65,100,153,156,157,158,159,160,200,201,202,285];
    
    # Create the object
    my $obj = Bio::Align::Subset->new(
                                      inputfile => $filename,
                                      format => $format
                                    );
    
    # View the result
    print Dumper($obj->build_subset($subset));

=cut

# Body ########################################################################
###############################################################################
###############################################################################
###############################################################################

=head1 CONSTRUCTOR

=head2 Bio::Align::Subset->new

=cut

###############################################################################
# Class data and methods
###############################################################################
{  
    # A list of all attributes wiht default values and read/write/required properties
    my %_attribute_properties = (
        _inputfile => ["????", "read.required"],
        _format    => ["????", "read.required"],
        _identifiers   => ["????", "read.write"   ],
        _sequences => ["????", "read.write"   ],
        _seq_length=> [0     , "read.write"   ]
    );
    
    # Global variable to keep count of existing objects
    my $_count = 0;
    
    # The list of all attributes
    sub _all_attributes {
        keys %_attribute_properties;
    }
    
    # Check if a given property is set for a given attribute
    sub _permissions{
        my ($self,$attribute, $permissions) = @_;
        $_attribute_properties{$attribute}[1] =~/$permissions/;
    }
    
    # Return the default value for a given attribute
    sub _attribute_default{
        my ($self,$attribute) = @_;
        $_attribute_properties{$attribute}[0];
    }
    
    # Manage the count of existing objects
    sub get_count{
        $_count;
    }
    sub _incr_count{
        ++$_count;
    }
    sub _decr_count{
        --$_count;
    }
    
}
#
# The constructor of the class
#
sub new {
    
    my ($class, %arg) = @_;
    my $self = bless {}, $class;
    
    foreach my $attribute ($self->_all_attributes()){
        
        # E.g. attribute = "_name", argument = "name"
        my ($argument) = ($attribute =~ /^_(.*)/);
        
        # If explicitly given
        if(exists $arg{$argument}){
            $self->{$attribute} = $arg{$argument};
        }
        
        # If not given but required
        elsif($self->_permissions($attribute, 'required')){
            croak("No $argument attribute as required");
        }
        
        # Set to default
        else{
            $self->{$attribute} = $self->_attribute_default($attribute);            
        }
        
    }
    
    # Called $class because it is a gobal method
    $class->_incr_count;
    
    $self->_extract_sequences;
    return $self;
    
}


#
# Obtaining the sequences in a Array
#
sub _extract_sequences{
    
    my $self = $_[0];
        
    my @identifiers;
    my @sequences;
    
    my $seqIO = Bio::SeqIO->new(
                             -file   => $self->get_inputfile,
                             -format => $self->get_format
                            );
    
    while( my $seq = $seqIO->next_seq){
        
        my $sequence_string = $seq->seq;
        $sequence_string =~ s/\s//g;
        
        push(@identifiers, $seq->id);
        $self->_verify_chain($sequence_string);
        push(@sequences, $sequence_string);
        
    }
    
    $self->set_identifiers(\@identifiers);
    $self->set_sequences(\@sequences);
    
}

#
# Build a subset
#
sub build_subset{
    
    my ($self, $subset) = @_;
    my %new_alignment = ();
    
    
    # Initialite array for the new sequences
    my @new_sequences = ();
    
    for(my $i=0;$i<=$#{$self->get_sequences};$i++){
        # Initialite a new string for the new sequence
        my $new_sequence = "";
        for my $index (@{$subset}){
            $new_sequence.= substr(${$self->get_sequences}[$i],$index,3);
        }
        push(@new_sequences, $new_sequence);
    }
    
    # Build the complete hash
    $new_alignment{sequences} = \@new_sequences;
    $new_alignment{headers}   = $self->get_identifiers;
    
    return \%new_alignment;
    
}

###############################################################################
# Auxiliary methods
###############################################################################
{
    #
    # Set the sequence length of the whole alignment
    #
    sub _set_sequence_length{
        my $self = $_[0];
        $self->{_seq_length} = $_[1];
    }
    
    #
    # Check if a the length of a given sequence match with the length of
    # the whole alignment.
    #
    sub _check_sequence_length{
        my $self = $_[0];
        my $tested_sequence_length = $_[1];
        $tested_sequence_length == $self->get_seq_length ? return 1 : return 0;
    }
    
    #
    # Verifies the integrity of a given sequence
    #
    sub _verify_chain{
        
        my ($self,$sequence) = @_;
        my $seq_length = length($sequence);
        
        
        # 1. The chain must be a DNA sequence
        $self->_isdna($sequence) ? 1 : carp("\nWARNING: The following sequence does not seems as a dna/rna (ATGCU) sequence:\n\n<< $sequence >>\n");
        
        # 2. Also, all the sequences must be equal. But if $_sequence_length
        # has not been updated, it takes the value of the length of this sequence.
        if($self->get_seq_length == 0){
            # The input file must be wrapped (non untermitated codons)
            $seq_length % 3 == 0 ? 1 : croak("The sequence length is not multiple of 3 ($seq_length)");
            $self->_set_sequence_length($seq_length);
        }else{
            $self->_check_sequence_length($seq_length) ? 1 : croak("A sequence length does not match with the length of the whole alignment");
        }
        return 1;
        
    }
    
    #
    # Verifies if a given string is a DNA sequence
    #
    sub _isdna{
        my ($self,$sequence) = ($_[0],uc($_[1]));
        if($sequence =~ /^[ACGTU]+$/){
             return 1;
        }else{
             return 0;
        }
    }
    
    
}
###############################################################################


###############################################################################
# Accessor Methods
###############################################################################
# This kind of method is called Accesor
# Method. It returns the value of a key
# and avoid the direct acces to the inner
# value of $obj->{_inputfile}.
###############################################################################
sub get_inputfile { $_[0] -> {_inputfile} }
sub get_format    { $_[0] -> {_format}    }
sub get_sequences { $_[0] -> {_sequences} }
sub get_identifiers   { $_[0] -> {_identifiers}   }
sub get_seq_length{ $_[0] -> {_seq_length}}
###############################################################################


###############################################################################
# Mutator Methods
###############################################################################
sub set_inputfile { my ($self, $inputfile) = @_;
                    $self-> {_inputfile} = $inputfile if $inputfile;
                  }
sub set_format    { my ($self, $format) = @_;
                    $self-> {_format} = $format if $format;
                  }
sub set_identifiers   { my ($self, $identifiers) = @_;
                    $self-> {_identifiers} = $identifiers if $identifiers;
                  }
sub set_sequences { my ($self, $sequences) = @_;
                    $self-> {_sequences} = $sequences if $sequences;
                  }
###############################################################################



# Footer ######################################################################
###############################################################################
###############################################################################
###############################################################################

=head1 AUTHOR

Hector Valverde and Juan Carlos Aledo, C<< <hvalverde@uma.es> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-align-subset at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Align-Subset>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Align::Subset


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Align-Subset>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Align-Subset>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Align-Subset>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Align-Subset/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2012 Hector Valverde and Juan Carlos Aledo.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of Bio::Align::Subset
