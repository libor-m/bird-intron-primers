@l = ();
push(@l, {'a' => 'b', 'c' => 'd'});
push(@l, {'a' => '1', 'c' => '2'});

#print "@l\n";
%h = $l[0];
print "val: " . $l[1]{'c'};