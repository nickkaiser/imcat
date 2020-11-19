# read arrays %ix{$chipno}, %ix{$chipno}, %cos{$chipno} from "config.db"


sub read_config_db {
	open (CONFIG_DB, 'config.db') || die "Can't open configuration database\n";
	while (<CONFIG_DB>) {
		unless (/^#/) {
			chop;
			@_ = split(' ');
	        	$chipno = shift @_;
			$ix{$chipno} = shift @_;
			$iy{$chipno} = shift @_;
			$cos{$chipno} = shift @_;
		}	
	}
}
1;


