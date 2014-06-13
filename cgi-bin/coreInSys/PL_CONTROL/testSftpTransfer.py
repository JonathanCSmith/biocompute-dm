#!/bin/env python


from sftpTransfer import sftpTransfer

ST=sftpTransfer()

#ST.username='brc00089'
#ST.loadSftpByUsername()
ST.seqProjectID=56
#ST.username='plavender'
ST.genPassword()
ST.genUsername()


ST.addSftpAccount()



#ST.delSftpAccount()

"""

ST.seqProjectID=89

"""
ST.setTransLocation()
ST.setTransfer()
ST.tranfer2sftp()

"""


#ST.delSftpAccount()


#ST.genPassword()
#ST.restoreSftpAccount()
#ST.addSftpAccount()



ST.seqProjectID=89
ST.genPassword()
ST.genUsername()

ST.addSftpAccount()

#ST.delSftpAccount()
"""



