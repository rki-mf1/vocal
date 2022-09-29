#!/usr/bin/env python3

# imports for mail functions
import argparse
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE
import os
from os.path import basename
import smtplib
import sys


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--msg-file", help="Text file with contents of message", default="-"
    )
    parser.add_argument(
        "--to", help="Email addresses to send message to (comma delimited)"
    )
    parser.add_argument(
        "--frm", help="Email addresses to send message from (one address)"
    )
    parser.add_argument("--cc", help="CC addresses (comma delimited)")
    parser.add_argument("--bcc", help="BCC addresses (comma delimited)")
    parser.add_argument("--subject", help="Subject")
    parser.add_argument("--attachments", help="Attachments", nargs="*")

    args = parser.parse_args()

    msg = ""
    if args.msg_file == "-":
        # Read message body from stdin
        msg = sys.stdin.read()
    else:
        if not os.path.exists(args.msg_file):
            print("Message file not found")
            return 1
        with open(args.msg_file) as f:
            msg = f.read()
        if args.msg_file.endswith(".html"):
            _html = True
        else:
            _html = False

    rkimail(
        msg, args.subject, args.to, args.frm, args.cc, args.bcc, args.attachments, _html
    )


# send internal rki mails
def rkimail(
    message,
    subject,
    to,
    frm="kongkitimanonk@rki.de",
    cc=None,
    bcc=None,
    files=None,
    html=False,
):
    mailserver = "smtp-connect.rki.local"
    mail = smtplib.SMTP(mailserver)
    # mail.ehlo()
    # max_limit_in_bytes = int( mail.esmtp_features['size'] )
    # sprint(max_limit_in_bytes)
    msg = MIMEMultipart()
    if html:
        msg.attach(MIMEText(message, "html"))
    else:
        msg.attach(MIMEText(message))
    msg["Subject"] = subject
    msg["From"] = frm

    if type(to) == str:
        msg["To"] = to
    elif type(to) == list:
        msg["To"] = COMMASPACE.join(to)
    else:
        raise TypeError("Wrong type. String or List are acceptable types!")

    if cc is not None:
        if type(cc) == str:
            msg["Cc"] = cc
        elif type(cc) == list:
            msg["Cc"] = COMMASPACE.join(cc)
        else:
            raise TypeError("Wrong type. String or List are acceptable types!")

    if bcc is not None:
        if type(bcc) == str:
            msg["Bcc"] = bcc
        elif type(bcc) == list:
            msg["Bcc"] = COMMASPACE.join(bcc)
        else:
            raise TypeError("Wrong type. String or List are acceptable types!")

    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(fil.read(), Name=basename(f))
        # After the file is closed
        part["Content-Disposition"] = 'attachment; filename="%s"' % basename(f)
        msg.attach(part)

    try:
        mail.send_message(msg)
    except Exception as err:
        raise
        return 1

    mail.quit()


# Main body
if __name__ == "__main__":
    main()
