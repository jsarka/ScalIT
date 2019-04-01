#!/bin/bash

binDir=/usr/home/wchen/osb/dosb/ndpr
wkDir=/usr/home/wchen/osb/dosb/ndpr

echo $binDir
echo $wkDir

$binDir/test_BR < $wkDir/Br.in > $wkDir/Br.out
$binDir/test_lr < $wkDir/lr.in > $wkDir/lr.out

