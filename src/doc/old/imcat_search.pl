#!/bin/perl

# this is my engin for searching the imcat html man pages
# please do not remove or change this fiel in any way.
# thanks, Nick

$webmaster = "Nick Kaiser (kaiser\@cita\.utoronto\.ca)";
$grep = "/bin/grep -i -l";
# $document_root = $ENV{'DOCUMENT_ROOT'};
$document_root = "../home/kaiser/imcatdoc";

&parse_form_data (*SEARCH);
$query = $SEARCH{'query'};

if ($query eq "") {
	&return_error (500, "Search Error", "Please enter a search query.");
} elsif ($query !~ /^(\w+)$/) {
	&return_error (500, "Search Error", "Invalid characters in query.");
} else {
	print "Content-type: text/html", "\n\n";
	print "<HTML>", "\n";
	print "<HEAD><TITLE>Imcat Man Pages Search Results</TITLE></HEAD>";
	print "<BODY>", "\n";
	print "<H1>Results of searching for imcat man pages for: ", $query, "</H1>";
	print "<HR>";

	$matches = 0;

	$file = "$document_root/mainindex.html";
#	&searchfile($file);

	open(LS, "/bin/ls -R $document_root |");
	while (<LS>) {
		chop($entry = $_);
		if (/^\.\./) {
			$entry =~ s/^$document_root\/(.*)/$1/;
			$entry =~ s/(.*):$/$1/;
			$dirname = $entry;
		} elsif (/\.html$/) {
			$basename = $entry;
			$fullname = "$document_root/$dirname/$basename";
			&searchfile($fullname);
		}
	}
	close(LS);

	
	print "<P>", "<HR>";
	print "Total number of matches: ", $matches, "<BR>";
	print "<HR>";
	print "</BODY></HTML>", "\n";

}

exit (0);

sub searchfile
{
	local ($source) = $_[0];
	local ($file);

#	print "searching $source <br>\n";
 
 	open (SEARCH, "$grep $query $source |");
	while (<SEARCH>) {
		chop($file = $_);
		$file =~ s/^$document_root\/(.*)/$1/;
		print qq|<A HREF="/~kaiser/imcatdoc/$file">$file</A><BR>|, "\n";
		$matches++;
	}
	close (SEARCH);
}


sub parse_form_data
{
    local (*FORM_DATA) = @_;

    local ( $request_method, $query_string, @key_value_pairs,
           $key_value, $key, $value);

    $request_method = $ENV{'REQUEST_METHOD'};

    if ($request_method eq "GET") {
        $query_string = $ENV{'QUERY_STRING'};
    } elsif ($request_method eq "POST") {
        read (STDIN, $query_string, $ENV{'CONTENT_LENGTH'});
    } else {
        &return_error (500, "Server Error",
                       "Server uses unsupported method");
    }

    @key_value_pairs = split (/&/, $query_string);

    foreach $key_value (@key_value_pairs) {
        ($key, $value) = split (/=/, $key_value);
        $value =~ tr/+/ /;
        $value =~ s/%([\dA-Fa-f][\dA-Fa-f])/pack ("C", hex ($1))/eg;

        if (defined($FORM_DATA{$key})) {
            $FORM_DATA{$key} = join ("\0", $FORM_DATA{$key}, $value);
        } else {
            $FORM_DATA{$key} = $value;
        }
    }
}

sub return_error
{
    local ($status, $keyword, $message) = @_;

    print "Content-type: text/html", "\n";
    print "Status: ", $status, " ", $keyword, "\n\n";

    print <<End_of_Error;

<title>CGI Program - Unexpected Error</title>
<h1>$keyword</h1>
<hr>$message</hr>
Please contact $webmaster for more information.

End_of_Error

    exit(1);
}
