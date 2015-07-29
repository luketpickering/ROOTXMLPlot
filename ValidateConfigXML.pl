#!/usr/bin/env perl

  use XML::SAX::ParserFactory;
  use XML::Validator::Schema;

 
  $num_args = $#ARGV + 1;
  if ($num_args != 1) {
    print "\nUsage: $0 <config.xml>\n";
    die;
  }
 
  $config_fn = $ARGV[0];
 
  print "Hello, $first_name $last_name\n";

  #
  # create a new validator object, using plots_config.xsd
  #
  $validator = XML::Validator::Schema->new(file => 'plots_config.xsd');

  #
  # create a SAX parser and assign the validator as a Handler
  #
  $parser = XML::SAX::ParserFactory->parser(Handler => $validator);

  #
  # validate foo.xml against foo.xsd
  #
  eval { $parser->parse_uri('$config_fn') };
  die "File failed validation: $@" if $@;
