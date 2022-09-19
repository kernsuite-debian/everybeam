#
# Regular cron jobs for the everybeam package
#
0 4	* * *	root	[ -x /usr/bin/everybeam_maintenance ] && /usr/bin/everybeam_maintenance
