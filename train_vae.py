import random
from pathlib import Path

import numpy as np
import matplotlib as plt
import torch
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
from torchvision.utils import make_grid
from tqdm import tqdm

import utils
import vae

# to ensure reproducible training/validation split
random.seed(42)

# find out if a GPU is available
if torch.cuda.is_available():
    device = torch.device("cuda")
elif torch.backends.mps.is_available():
    device = torch.device("mps")
else:
    device = torch.device("cpu")

# directorys with data and to store training checkpoints and logs
DATA_DIR = Path.cwd().parent / "TrainingData"
CHECKPOINTS_DIR = Path.cwd() / "vae_model_weights"
CHECKPOINTS_DIR.mkdir(parents=True, exist_ok=True)
TENSORBOARD_LOGDIR = "vae_runs"

# training settings and hyperparameters
NO_VALIDATION_PATIENTS = 2
IMAGE_SIZE = [64, 64]
BATCH_SIZE = 32
N_EPOCHS = 200
DECAY_LR_AFTER = 50
LEARNING_RATE = 1e-4
DISPLAY_FREQ = 10

# dimension of VAE latent space
Z_DIM = 256

# function to reduce the
def lr_lambda(the_epoch):
    """Function for scheduling learning rate"""
    return (
        1.0
        if the_epoch < DECAY_LR_AFTER
        else 1 - float(the_epoch - DECAY_LR_AFTER) / (N_EPOCHS - DECAY_LR_AFTER)
    )


# find patient folders in training directory
# excluding hidden folders (start with .)
patients = [
    path
    for path in DATA_DIR.glob("*")
    if not any(part.startswith(".") for part in path.parts)
]
random.shuffle(patients)

# split in training/validation after shuffling
partition = {
    "train": patients[:-NO_VALIDATION_PATIENTS],
    "validation": patients[-NO_VALIDATION_PATIENTS:],
}

# load training data and create DataLoader with batching and shuffling
dataset = utils.ProstateMRDataset(partition["train"], IMAGE_SIZE)
dataloader = DataLoader(
    dataset,
    batch_size=BATCH_SIZE,
    shuffle=True,
    drop_last=True,
    pin_memory=True,
)

# load validation data
valid_dataset = utils.ProstateMRDataset(partition["validation"], IMAGE_SIZE)
valid_dataloader = DataLoader(
    valid_dataset,
    batch_size=BATCH_SIZE,
    shuffle=True,
    drop_last=True,
    pin_memory=True,
)

# initialise model, optimiser
vae_model = vae.VAE(IMAGE_SIZE, Z_DIM).to(device)
optimizer = torch.optim.Adam(vae_model.parameters(), lr=LEARNING_RATE)
# add a learning rate scheduler based on the lr_lambda function
scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lr_lambda)

# training loop
writer = SummaryWriter(log_dir=TENSORBOARD_LOGDIR)  # tensorboard summary
for epoch in range(N_EPOCHS):
    current_train_loss = 0.0
    current_valid_loss = 0.0
    
    # training iterations
    # the required implementation of training iterations in pytorch is composed of 5 steps:
    for inputs, labels in tqdm(dataloader, position=0):
        optimizer.zero_grad() # (1) zeroing the gradients in each iterations
        outputs, mu, logvar = vae_model(inputs.to(device))  # (2) forward pass of the model
        loss = vae.vae_loss(inputs.to(device), outputs, mu, logvar) # (3) computing the loss
        loss.backward()  # (4) backpropagating the loss
        current_train_loss += loss.item()
        optimizer.step()  # (5) stepping the optimiser (update the weights)


    # evaluate validation loss
    with torch.no_grad():
        vae_model.eval()
        for inputs, labels in tqdm(valid_dataloader, position=0):
            outputs, mu, logvar = vae_model(inputs.to(device))  # (2) forward pass of the model
            loss = vae.vae_loss(inputs.to(device), outputs, mu, logvar) # (3) computing the loss
            current_valid_loss += loss.item()
            
        vae_model.train()

    # write to tensorboard log
    writer.add_scalar("Loss/train", current_train_loss / len(dataloader), epoch)
    writer.add_scalar(
        "Loss/validation", current_valid_loss / len(valid_dataloader), epoch
    )
    scheduler.step() # step the learning step scheduler

    # save examples of real/fake images
    x_recon = outputs.cpu()
    x_real = inputs.cpu()
    if (epoch + 1) % DISPLAY_FREQ == 0:
        img_grid = make_grid(
            torch.cat((x_recon[:5], x_real[:5])), nrow=5, padding=12, pad_value=-1
        )
        writer.add_image(
            "Real_fake", np.clip(img_grid[0][np.newaxis], -1, 1) / 2 + 0.5, epoch + 1
        )
    
    noise = vae.get_noise(BATCH_SIZE, Z_DIM, device)# sample noise 
    with torch.no_grad():
         vae_model.eval()
         generated_images = vae_model.generator(noise) # generate images and display
         vae_model.train()
    generated_image_grid = make_grid(generated_images, nrow=8, padding=2, pad_value=1)
    writer.add_image("GeneratedImages", generated_image_grid, epoch+1)

torch.save(vae_model.state_dict(), CHECKPOINTS_DIR / "vae_model.pth")