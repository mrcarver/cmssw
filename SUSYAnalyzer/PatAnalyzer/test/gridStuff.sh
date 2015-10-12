voms-proxy-init -cert ~/.globus/usercert.pem -key ~/.globus/userkey.pem -out ~/user.proxy -voms cms -valid 168:30
voms-proxy-init
export X509_USER_PROXY=$(voms-proxy-info -path 2>/dev/null)
cp /tmp/x509up_u5278 testGrid
