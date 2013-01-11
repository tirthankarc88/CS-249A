$mypath = `pwd`; 
chomp($mypath);

`make clean && make`;

$testpath = "../test-cases/";
#$testpath = $mypath;

chdir($testpath) or die "$!";

@testfiles = `ls *.log`;
chomp(@testfiles);

chdir($mypath) or die "$!";

foreach $file(@testfiles){
    print "Running $file ... ";
    $result = `./asgn1 $testpath/$file`;
    $ref = `cat $testpath/$file.ref2`;
    #$ref = `cat ../$file.out`;
    if ($ref eq $result) {
	print "Test passed!\n";
    } else {
	print "Test failed!\n";
    }
}
