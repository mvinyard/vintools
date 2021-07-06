def _initialize_gmail_server(sender, email_password):
    
    """"""
    import smtplib
    # creates SMTP session
    server = smtplib.SMTP( "smtp.gmail.com", 587 ) # 587 for gmail
    # start TLS for security
    server.starttls()
    server.login(sender, email_password)
    
    return server

def _prepare_message(content, recipient, subject, cc, sender):
    
    from email.message import EmailMessage

    msg = EmailMessage()   
    msg["Subject"] = subject
    msg["From"] = "michaele.vinyard@gmail.com"
    msg["To"] = recipient
    msg["Cc"] = cc
    
    msg.set_content(content)
        
    return msg
    
def _send_email(content, email_password, subject=None, recipient=None, cc=None, sender='michaele.vinyard@gmail.com'):
    
    """
    Code to send a simple email. 
    
    Parameters:
    -----------
    content
        Text email contents.
    
    email_password
        password to email server being used.
    
    subject
        email subject line to be included.
        [optional] default: None
    
    recipient
        email address to recieve the message.
        [optional] default: None
    
    cc
        email address to recieve the message. carbon-copied.
        [optional] default: None
    
    sender
        email address from which emails should be sent. 
        Must match server being used (gmail) and the password provided.
        [optional] default: michaele.vinyard@gmail.com
    
    Returns:
    --------
    None
        email is sent. 
    """
    
    server = _initialize_gmail_server(sender, email_password)
    msg = _prepare_message(content, recipient, subject, cc, sender)

    server.send_message(msg)
    server.quit()
    
    print("email sent to:", recipient, cc)
