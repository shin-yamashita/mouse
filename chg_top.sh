

for fn in ls html/*.html ; do
 dst=`echo $fn|sed 's;^html/;;g'`
 echo $fn, $dst
 cat $fn | sed 's|http://localhost/~shin|.|g;s|../cache/|image/|g;s|../image/|image/|g;s|../refimg/|refimg/|g;s|../skin/|skin/|g' > $dst
done

